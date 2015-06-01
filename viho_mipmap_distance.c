// Homography viewer.
//
// This program loads a color image and allows the user to drag the four
// corners of the image around the window, thus deforming the image by an
// homography.
//
//
// Usage:
//
// 	viho image.png
//
//
// Mouse controls:
//
// 	1. Drag the control points with the left mouse button to change their
// 	position in the window
//
//	2. Drag the control points with the right mouse button to change their
//	position in the image
//
//	3. Mouse wheel for zoom-in and zoom-out
//
//
// Keys:
//
//	q	exit the viewer
//	c	reset viewer and center image
//
// 	+	zoom-in
// 	-	zoom-out
// 	ARROWS	move the view
//
// 	0	use nearest-neighbor interpolation
// 	1	use linear interpolation
// 	2	use bilinear interpolation
// 	3	use bicubic interpolation
//
// 	p	toggle periodic extrapolation
// 	w	toggle horizon
// 	.	toggle grid points
//
//Distance
// default: plus grand cote du prgm 
// 1: prgm plus grande diag
// 2: prgm plus petite diag
// 3: valeur singulière plus grande
// 4: valeur singulière plus petite
// 5: determinant de l'affinité



// Compilation:
//
//	c99 -O3 -DNDEBUG viho.c iio.c ftr.c -lX11 -lpng -ljpeg -ltiff -o viho 
//  ou plutot :
//  clang -I/usr/X11/include viho.c iio.c ftr.c -L/usr/X11/lib -lX11


// SECTION 1. Libraries and data structures                                 {{{1

// standard libraries
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

// user interface library
#include "ftr.h"


// radius of the disks that are displayed around control points
#define DISK_RADIUS 7

// zoom factor for zoom-in and zoom-out
#define ZOOM_FACTOR 1.43


// data structure to store the state of the viewer
struct viewer_state {
	// image data
	float *img;
	int iw, ih, pd;

	// geometry
	double c[4][2]; // control points in window coordinates
	double p[4][2]; // control points in image coordinates
	double offset[2], scale; // window viewport

	// dragging state
	bool dragging_point;
	bool dragging_ipoint;
	int dragged_point;

	// display options
	int interpolation_order; // 0=nearest, 1=linear, 2=bilinear, 3=bicubic
	bool tile_plane;
	bool show_horizon;
	bool show_grid_points;
};


// function to reset and center the viewer
static void center_view(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// set control points near the corners of the window
	double margin = 33;
	e->c[0][0] = margin;
	e->c[0][1] = margin;
	e->c[1][0] = f->w - margin;
	e->c[1][1] = margin;
	e->c[2][0] = f->w - margin;
	e->c[2][1] = f->h - margin;
	e->c[3][0] = margin;
	e->c[3][1] = f->h - margin;

	// set control points exactly at the corners of the image
	e->p[0][0] = 0;
	e->p[0][1] = 0;
	e->p[1][0] = e->iw - 1;
	e->p[1][1] = 0;
	e->p[2][0] = e->iw - 1;
	e->p[2][1] = e->ih - 1;
	e->p[3][0] = 0;
	e->p[3][1] = e->ih - 1;

	// drag state
	e->dragging_point = false;
	e->dragging_ipoint = false;

	// viewport
	e->offset[0] = 0;
	e->offset[1] = 0;
	e->scale = 1;

	// visualization options
	e->interpolation_order = 0;
	e->tile_plane = false;
	e->show_horizon = false;
	e->show_grid_points = false;
}


// funtion to test whether a point is inside the window
static int insideP(struct FTR *f, int x, int y)
{
	return x >= 0 && y >= 0 && x < f->w && y < f->h;
}



// SECTION 2. Linear algebra                                                {{{1

// y = H(x)
static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}

// compute the inverse homography (inverse of a 3x3 matrix)
static double invert_homography(double invH[3][3], double H[3][3])
{
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = (a[4]*a[8]-a[5]*a[7])/det;
	r[1] = (a[2]*a[7]-a[1]*a[8])/det;
	r[2] = (a[1]*a[5]-a[2]*a[4])/det;
	r[3] = (a[5]*a[6]-a[3]*a[8])/det;
	r[4] = (a[0]*a[8]-a[2]*a[6])/det;
	r[5] = (a[2]*a[3]-a[0]*a[5])/det;
	r[6] = (a[3]*a[7]-a[4]*a[6])/det;
	r[7] = (a[1]*a[6]-a[0]*a[7])/det;
	r[8] = (a[0]*a[4]-a[1]*a[3])/det;
	return det;
}

// C = AoB, composition of two homographies (product of 3x3 matrices)
static void compose_homographies(double C[3][3], double A[3][3], double B[3][3])
{
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		C[i][j] = 0;
		for (int k = 0; k < 3; k++)
			C[i][j] += A[i][k] * B[k][j];
	}
}

// Find the homography that changes the canonical projective basis
// into the given four points (x, y, z, w)
static void homography_from_four_points(double H[3][3],
		double x[2], double w[2], double z[2], double y[2])
{
	// fix the degree of freedom (assuming the four points are finite)
	double t = 1;

	// translation coefficients
	double p = x[0];
	double q = x[1];

	// "core" 2x2 system
	double A = w[0] - z[0];
	double B = y[0] - z[0];
	double C = w[1] - z[1];
	double D = y[1] - z[1];
	double P = z[0] - y[0] - w[0] + p;
	double Q = z[1] - y[1] - w[1] + q;
	double DET = A * D - B * C;
	double r = (D * P - B * Q) / DET;
	double s = (A * Q - C * P) / DET;
	if (!isnormal(DET))
		fprintf(stderr, "denormal! DET = %g\n", DET);

	// solve the rest of the diagonal system
	double a = w[0] * ( 1 + r ) - p;
	double b = y[0] * ( 1 + s ) - p;
	double c = w[1] * ( 1 + r ) - q;
	double d = y[1] * ( 1 + s ) - q;

	// fill-in the output
	H[0][0] = a; H[0][1] = b; H[0][2] = p;
	H[1][0] = c; H[1][1] = d; H[1][2] = q;
	H[2][0] = r; H[2][1] = s; H[2][2] = t;
}

// Find the homography that moves the four points (x,y,z,w) to (a,b,c,d)
static void homography_from_eight_points(double H[3][3],
		double x[2], double y[2], double z[2], double w[2],
		double a[2], double b[2], double c[2], double d[2])
{
	double H1[3][3], H2[3][3], iH1[3][3];
	homography_from_four_points(H1, x, y, z, w);
	homography_from_four_points(H2, a, b, c, d);
	invert_homography(iH1, H1);
	compose_homographies(H, H2, iH1);
}

/// compute the vector product of two vectors
static void vector_product(double axb[3], double a[3], double b[3])
{
	// a0 a1 a2
	// b0 b1 b2
	axb[0] = a[1] * b[2] - a[2] * b[1];
	axb[1] = a[2] * b[0] - a[0] * b[2];
	axb[2] = a[0] * b[1] - a[1] * b[0];
}




// SECTION 3. Coordinate Conversions                                        {{{1

// This program deals with three systems of coordinates
//
// "image"  : coordinates in the rectangular domain of the input image
// "view"   : coordinates in the infinite plane where the image is mapped
// "window" : coordinates in the window, which is a rectangluar piece of "view"
//

// change from view coordinates to window coordinates
static void map_view_to_window(struct viewer_state *e, double y[2], double x[2])
{
	for (int k = 0; k < 2; k++)
		y[k] = e->offset[k] + e->scale * x[k];
}

// change from window coordinates to view coordinates
static void map_window_to_view(struct viewer_state *e, double y[2], double x[2])
{
	for (int k = 0; k < 2; k++)
		y[k] = ( x[k] - e->offset[k] ) / e->scale;
}

// obtain the direct homography from the current configuration
static void obtain_current_homography(double H[3][3], struct viewer_state *e)
{
	double C[4][2];
	for (int p = 0; p < 4; p++)
		map_view_to_window(e, C[p], e->c[p]);
	homography_from_eight_points(H,
			C[0], C[1], C[2], C[3],
			e->p[0], e->p[1], e->p[2], e->p[3]
			);
}

// change from window coordinates to image coordinates
static void map_window_to_image(struct viewer_state *e, double *y, double *x)
{
	double H[3][3], C[4][2];
	obtain_current_homography(H, e);
	apply_homography(y, H, x);
}

// Convert a floating-point color into a byte in the range [0,255]
// (Colors outisde of this range are saturated.)
static double float_to_byte(double x)
{
	if (x < 0) return 0;
	if (x < 255) return x;
	return 255;
}


// SECTION 4. Extrapolation                                                 {{{1

// A "extrapolator" evaluates an image at an arbitrary integral position.
// When the position is outside the image domain, the value is extrapolated
// (by periodization, or by a constant value).

// type of the "extrapolator" functions
typedef float (*extrapolator_t)(float*,int,int,int,int,int,int);

// auxiliary function: compute n%p correctly, even for huge and negative numbers
static int good_modulus(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	return r;
}

// instance of "extrapolator_t", extrapolate by periodicity
static float getsample_per(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// instance of "extrapolator_t", extrapolate by a constant value
static float getsample_cons(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return value;
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// obtain an extrapolator function from the current configuration
static extrapolator_t obtain_extrapolator(struct viewer_state *e)
{
	if (e->tile_plane)
		return getsample_per;
	float top = 255;
	getsample_cons(&top, 0, 0, 0, 0, 0, 0);
	return getsample_cons;
}



// SECTION 5. Local Interpolation                                           {{{1
//
// An "interpolator" evaluates an image at an arbitrary floating-point
// position.  When the position is not integral, the value is a combination of
// neighboring values.  Notice that when the position is outside the image
// domain, an extrapolator is also needed.

// type of the "interpolator" functions
typedef float (*interpolator_t)(float*,int,int,int,float,float,int,
		extrapolator_t);

// auxiliary function for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}

// instance of "interpolator_t", for bilinear interpolation
static float bilinear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = pix(x, w, h, pd, ip  , iq  , l);
	float b = pix(x, w, h, pd, ip+1, iq  , l);
	float c = pix(x, w, h, pd, ip  , iq+1, l);
	float d = pix(x, w, h, pd, ip+1, iq+1, l);
	return evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
}

// auxiliary code for "linear" interpolation
#include "marchi.c"

// instance of "interpolator_t", for linear (marching) interpolation
static float linear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = pix(x, w, h, pd, ip  , iq  , l);
	float b = pix(x, w, h, pd, ip+1, iq  , l);
	float c = pix(x, w, h, pd, ip  , iq+1, l);
	float d = pix(x, w, h, pd, ip+1, iq+1, l);
	return marchi(a, c, b, d, p-ip, q-iq);
}

// instance of "interpolator_t" for nearest neighbor interpolation
static float nearest_neighbor_at(float *x, int w, int h, int pd,
		float p, float q, int l, extrapolator_t pix)
{
	int ip = round(p);
	int iq = round(q);
	return pix(x, w, h, pd, ip, iq, l);
}


// one-dimensional cubic interpolation of four data points ("Keys")
static float cubic_interpolation(float v[4], float x)
{
	return v[1] + 0.5 * x*(v[2] - v[0]
			+ x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
			+ x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

// two-dimensional separable cubic interpolation, on a 4x4 grid
static float bicubic_interpolation_cell(float p[4][4], float x, float y)
{
	float v[4];
	v[0] = cubic_interpolation(p[0], y);
	v[1] = cubic_interpolation(p[1], y);
	v[2] = cubic_interpolation(p[2], y);
	v[3] = cubic_interpolation(p[3], y);
	return cubic_interpolation(v, x);
}

// instance of "interpolator_t" for bicubic interpolation
static float bicubic_interpolation_at(float *img, int w, int h, int pd,
		float x, float y, int l, extrapolator_t p)
{
	x -= 1;
	y -= 1;

	int ix = floor(x);
	int iy = floor(y);
	float c[4][4];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			c[i][j] = p(img, w, h, pd, ix + i, iy + j, l);
	return bicubic_interpolation_cell(c, x - ix, y - iy);
}


// obtain an interpolation operator from the current configuration
static interpolator_t obtain_interpolator(struct viewer_state *e)
{
	if (e->dragging_point || e->dragging_ipoint) return nearest_neighbor_at;
	if (e->interpolation_order == 1) return linear_interpolation_at;
	if (e->interpolation_order == 2) return bilinear_interpolation_at;
	if (e->interpolation_order == 3) return bicubic_interpolation_at;
	return nearest_neighbor_at;
}



//Mip-Map section by Momo

int coo (int i,int j,int w, int pd){
	return (i+j*w)*pd;
}


//on périodise par défaut
int coord_mipmap (int u,int v,int w,int pd,int d,int l){
	int x,y;
	int z = w/pow(2,d-1);
	x=good_modulus(u,z);
	y=good_modulus(v,z);
	return 2*w*w*(1-1/pow(2,d-1))*pd+(x+y*w)*pd+l;
}

//pour plus d'option, on déclare ce filtre qui est la moyenne
//float filter_mipmap(int i,int j,int w,int pd,int d,int l,float *r){
//	return (r[coord_mipmap(2*i,2*j,w,pd,d-1,l)]+r[coord_mipmap(2*i+1,2*j,w,pd,d-1,l)]+r[coord_mipmap(2*i,2*j+1,w,pd,d-1,l)]+r[coord_mipmap(2*i+1,2*j+1,w,pd,d-1,l)])/4;
//}

#define SIGC 0.36 //paramètre du filtre gaussien au carre
#define TAPS 5 //racine du nombre de coefficient interpolé, doit être impaire

float filter_mipmap_aux(int u,int v,int w,int pd,int d,int l,float *r,double *g){
	float total=0;
	int i,j;
	for (i=0;i<TAPS;i++){
		for(j=0;j<TAPS;j++){
			total = total + g[i+TAPS*j]*r[coord_mipmap(u-(TAPS-1)/2+i,v-(TAPS-1)/2+j,w,pd,d-1,l)];
		}
	}
	return total;
}

float filter_mipmap(int u,int v,int w,int pd,int d,int l,float *r,double *g){
	float a[4];
	a[0] = filter_mipmap_aux(u,v,w,pd,d,l,r,g);
	a[1] = filter_mipmap_aux(u+1,v,w,pd,d,l,r,g);
	a[2] = filter_mipmap_aux(u,v+1,w,pd,d,l,r,g);
	a[3] = filter_mipmap_aux(u+1,v+1,w,pd,d,l,r,g);
	return (a[0]+a[1]+a[2]+a[3])/4;
}

// on admet pour l'instant que l'image et de taille une puissance de 2, cette fonction construit le mip_map d'une image
float *mip_map(float *img,int w,int h, int pd, float *r){
int i,j,d,l,ll;
double gauss[TAPS*TAPS];
for (i=0;i<TAPS;i++){
	for(j=0;j<TAPS;j++){
		gauss[i+TAPS*j]=exp(-(pow(j-(TAPS-1)/2,2)+pow(i-(TAPS-1)/2,2))/(2*SIGC))/(2*M_PI*SIGC);   //on a prit sigma = 0,6 (ptet le prendre variable)
	}
}
float norm = 0;
for(i=0;i<pow(TAPS,2);i++){norm = norm + gauss[i];}
for(i=0;i<pow(TAPS,2);i++){gauss[i]=gauss[i]/norm;}
for(ll=0;ll<3;ll++){
	if(ll>pd-1){l=pd-1;}else{l=ll;}
	for(i=0;i<w;i++){
		for(j=0;j<h;j++){
			r[(i+j*w)*3+ll]=img[(i+j*w)*pd+l];
		}
	}
}
for(l=0;l<3;l++){
	for(d=2;pow(2,d-1)<=w;d++){
		for(j=0;j<w/pow(2,d-1);j++){
			for(i=0;i<w/pow(2,d-1);i++){
				r[coord_mipmap(i,j,w,3,d,l)]= filter_mipmap(2*i,2*j,w,3,d,l,r,gauss);			
			}
		}
	}
}
return r;
}



// on a utilisé cela comme ref pour écrire les fonctions
// du/dx = D[1]y+D[2]
// du/dy = D[5]x+D[6]
// dv/dx = D[3]y+D[4]
// dv/dy = D[7]x+D[8]

void precal_D_vsing(double *D){
	//et pour E le coeff dans le polynôme
	D[12]=D[5]*D[7];
	D[13]=D[6]*D[7]+D[5]*D[8];
	D[14]=D[1]*D[3];
	D[15]=D[1]*D[4]+D[3]*D[2];
	D[16]=D[6]*D[8]+D[2]*D[4];
}


//D(x,y)=sqrt(12xy+13x+14y+15)/c^2
void precal_D_det(double *D){
	D[12]=D[1]*D[7]+D[5]*D[3];
	D[13]=D[2]*D[7]+D[5]*D[4];
	D[14]=D[3]*D[6]+D[1]*D[8];
	D[15]=D[2]*D[8]+D[6]*D[4];
}


void precal_D_sqr(double H[3][3],double *D){
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0];
	D[1]=a[0]*a[7]-a[6]*a[1];
	D[2]=a[0]*a[8]-a[6]*a[2];
	D[3]=a[3]*a[7]-a[6]*a[4];
	D[4]=a[3]*a[8]-a[6]*a[5];
	D[5]=-D[1];
	D[6]=a[1]*a[8]-a[7]*a[2];
	D[7]=-D[3];
	D[8]=a[4]*a[8]-a[7]*a[5];
	D[9]=a[6];
	D[10]=a[7];
	D[11]=a[8];
}

void precal_D_prgm(double H[3][3],double *D){
	double *a = H[0];
	D[1]=a[0]*a[7]-a[6]*a[1]-(a[3]*a[7]-a[6]*a[4]);
	D[2]=a[0]*a[8]-a[6]*a[2]-(a[3]*a[8]-a[6]*a[5]);
	D[3]=a[3]*a[7]-a[6]*a[4]-(a[0]*a[7]-a[6]*a[1]);
	D[4]=a[1]*a[8]-a[7]*a[2]-(a[4]*a[8]-a[7]*a[5]);
	D[5]=a[0]*a[7]-a[6]*a[1]+(a[3]*a[7]-a[6]*a[4]);
	D[6]=a[0]*a[8]-a[6]*a[2]+(a[3]*a[8]-a[6]*a[5]);
	D[7]=a[3]*a[7]-a[6]*a[4]+(a[0]*a[7]-a[6]*a[1]);
	D[8]=a[1]*a[8]-a[7]*a[2]+(a[4]*a[8]-a[7]*a[5]);	
	D[9]=a[6];
	D[10]=a[7];
	D[11]=a[8];
}


int precal_D(double H[3][3],double *D,int fun_dist){
	if(fun_dist==1 || fun_dist==2){precal_D_prgm(H,D); return 0;}
	if(fun_dist==3 || fun_dist==4){precal_D_sqr(H,D); precal_D_vsing(D); return 0;}
	if(fun_dist==5){precal_D_sqr(H,D); precal_D_det(D); return 0;}
	precal_D_sqr(H,D); return 0;
}



// on note D = (1,2,3,4,5,6,7,8,9,10,11) avec 
// D = max(sqrt( (1y+2)^2 + (3y+4)^2 )),sqrt( (5x+6)^2 + (7x+8)^2 ))/(9x + 10y +11)^2

//lien 1 = 1-3
//2 = 2-4
//3 = 5-7
//4 = 6-8

//et de même avec plus

//Distance variables par defaut max_sqr
//1: prgm diag    max ( sqrt ( (1y+2)^2 + (3x+4)^2 ) , sqrt ( (5y+6)^2 + (7x+8)^2 )/(9x + 10y +11)^2

#define D_BIAS 0.

double absd(double p){if(p>0){return p;}{return -p;}}

double cal_D_det(int x,int y,double *D){
	double p[2]={x,y};
	double a,b,c;
	c = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	a = sqrt(absd(D[12]*p[0]*p[1]+D[13]*p[0]+D[14]*p[1]+D[15]));
	return(a/c+D_BIAS);	
}

// elle prend des int est retourne un double
double cal_D_sqr(int x,int y,double *D){
	double p[2]={x,y};
	double a,b,c;
	c = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	c = c;
	a = sqrt (pow(D[1]*p[1]+D[2],2) + pow(D[3]*p[1]+D[4],2));
	b = sqrt (pow(D[5]*p[0]+D[6],2) + pow(D[7]*p[0]+D[8],2));
	if(a>b){return a/c+D_BIAS;}{return b/c+D_BIAS;}
}

//diff subtil, c'est du code pour rien
double cal_D_prgm(int x,int y,double *D,int opt){
	double p[2]={x,y};
	double a,b,c;
	c = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	c = c*sqrt(2);
	a = sqrt (pow(D[1]*p[1]+D[2],2) + pow(D[3]*p[0]+D[4],2));
	b = sqrt (pow(D[5]*p[1]+D[6],2) + pow(D[7]*p[0]+D[8],2));
	if(opt==1){if(a>b){return a/c+D_BIAS;}{return b/c+D_BIAS;}}{if(a<b){return a/c+D_BIAS;}{return b/c+D_BIAS;}}
}

double cal_D_vsing(int x,int y,double *D,int opt){
	double p[2]={x,y};
	double A,B,C,E;
	C = pow(D[9]*p[0]+D[10]*p[1]+D[11],2)*sqrt(2);  //(car le /2 dans le pol de deg 2)
	A = pow(D[1]*p[1]+D[2],2) + pow(D[3]*p[1]+D[4],2);
	B = pow(D[5]*p[0]+D[6],2) + pow(D[7]*p[0]+D[8],2);
	E = pow(D[12]*pow(p[0],2)+D[13]*p[0]+D[14]*pow(p[1],2)+D[15]*p[1]+D[16],2);    // E = du/dx dv/dx + du/dy dv/dy, E = (12x^2+13x+14y^2+15y+16)^2
	if(opt==3){
		return sqrt(A+B+sqrt(pow(A+B,2)-4*(A*B-E)))/C + D_BIAS;
	}{
		return sqrt(A+B-sqrt(pow(A+B,2)-4*(A*B-E)))/C + D_BIAS;
	} //et le moins pour avoir lambda au lieu de t //sqrt(A+B-sqrt(pow(A+B,2)-4*(AB-E)))/(2*c)
}


double cal_D(int x,int y,double *D,int fun_dist){
	if(fun_dist==1 || fun_dist==2){return cal_D_prgm(x,y,D,fun_dist);}
	if(fun_dist==3 || fun_dist==4){return cal_D_vsing(x,y,D,fun_dist);}
	if(fun_dist==5){return cal_D_det(x,y,D);}
	return cal_D_sqr(x,y,D);
}


//int coord_mipmap (int u,int v,int w, int pd,int d,int l)
// for bilinear interpolation during mipmap
static float bilinear_mipmap(float *x, int w, int h,
		float p, float q, int l, extrapolator_t pix, int d)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = x[coord_mipmap(ip,iq,w,3,d,l)];
	float b = x[coord_mipmap(ip+1,iq,w,3,d,l)];
	float c = x[coord_mipmap(ip,iq+1,w,3,d,l)];
	float dd = x[coord_mipmap(ip+1,iq+1,w,3,d,l)];
	return evaluate_bilinear_cell(a, b, c, dd, p-ip, q-iq);
}



//mipmapping interpolator, img is a mipmapping
static float mip_mapping_interpolator_at(float *img, int w, int h, int pd,
		float x, float y, int l, extrapolator_t p,double D,int opt){
	int d ;
	float a,b;
	for(d=0;(pow(2,d))<=D && (pow(2,d))<=w;d++){;}
	if (d==1||d==0){return bilinear_mipmap(img,w,h,x,y,l,p,1);}
	if ((pow(2,d))>w){return img[6*w*(w-1)+l];} 
	a = bilinear_mipmap(img,w,h,x/pow(2,d-1),y/pow(2,d-1),l,p,d);
	b = bilinear_mipmap(img,w,h,x/pow(2,d-2),y/pow(2,d-2),l,p,d-1);
//	if(opt == 7){return ((D-pow(2,d-1))/(pow(2,d-1)))*a + ((pow(2,d)-D)/pow(2,d-1))*b;}
	return (log2(D)-d+1)*a + (d-log2(D))*b;  //((D-pow(2,d-1))/(pow(2,d-1)))*a + ((pow(2,d)-D)/pow(2,d-1))*b; //logarithmique ? (log(D)-d+1)*a + (d-log(D))*b;
}


// on a l'impression que img[i][j][l] avec l la couleur, pd le nombre de couleur = img[(i+j*w)*pd + l]
// et dans le mip-map ?
// r[u][v][d] = h*(w+..+w/(2^(d-1)))*pd + (u+v*w/d)*pd+l


// SECTION 6. Main Warping Function                                         {{{1

// draw the image warped by the current homography
static void draw_warped_image(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	int w = e->iw;
	int h = e->ih;

	double         H[3][3];   obtain_current_homography(H, e);

//la fonction de distance correspondante à l'homographie (double)
	double D[17];
	int fun_dist = e->interpolation_order;
	precal_D(H,D,fun_dist);
	extrapolator_t OUT      = obtain_extrapolator(e);
//	interpolator_t EVAL     = obtain_interpolator(e); //on périodise par défaut

	for (int j = 0; j < f->h; j++)
	for (int i = 0; i < f->w; i++)
	{
		double p[2] = {i, j};
		double d = cal_D(i,j,D,fun_dist);
		apply_homography(p, H, p);
		p[0] = (p[0] - 0.5) * w / (w - 1.0);
		p[1] = (p[1] - 0.5) * h / (h - 1.0);
		for (int l = 0; l < 3; l++)
		{
			int idx = l + 3 * (f->w * j + i);
			float v = mip_mapping_interpolator_at(e->img, w, h, e->pd, p[0], p[1], l, OUT, d, fun_dist);
			f->rgb[idx] = float_to_byte(v);
		}
	}
}



// SECTION 7. Drawing                                                       {{{1

// Subsection 7.1. Drawing segments                                         {{{2

// generic function to traverse a segment between two pixels
void traverse_segment(int px, int py, int qx, int qy,
		void (*f)(int,int,void*), void *e)
{
	if (px == qx && py == qy)
		f(px, py, e);
	else if (qx + qy < px + py) // bad quadrants
		traverse_segment(qx, qy, px, py, f, e);
	else {
		if (qx - px > qy - py || px - qx > qy - py) { // horizontal
			float slope = (qy - py)/(float)(qx - px);
			for (int i = 0; i < qx-px; i++)
				f(i+px, lrint(py + i*slope), e);
		} else { // vertical
			float slope = (qx - px)/(float)(qy - py);
			for (int j = 0; j <= qy-py; j++)
				f(lrint(px+j*slope), j+py, e);
		}
	}
}

// auxiliary function for drawing a red pixel
static void plot_pixel_red(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 255;
		f->rgb[3*idx+1] = 0;
		f->rgb[3*idx+2] = 0;
	}
}

// auxiliary function for drawing a green pixel
static void plot_pixel_green(int x, int y, void *e)
{
	struct FTR *f = e;
	if (insideP(f, x, y)) {
		int idx = f->w * y + x;
		f->rgb[3*idx+0] = 0;
		f->rgb[3*idx+1] = 128;
		f->rgb[3*idx+2] = 0;
	}
}

// function to draw a red segment
static void plot_segment_red(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_red, f);
}

// function to draw a green segment
static void plot_segment_green(struct FTR *f,
		double x0, double y0, double xf, double yf)
{
	traverse_segment(x0, y0, xf, yf, plot_pixel_green, f);
}



// Subsection 7.2. Drawing user-interface elements                          {{{2

// draw the positions of the image samples
static void draw_grid_points(struct FTR *f)
{
	struct viewer_state *e = f->userdata;
	double H[3][3], iH[3][3];
	obtain_current_homography(iH, e);
	invert_homography(H, iH);
	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 0; // black
	for (int j = 0; j < e->ih; j++)
	for (int i = 0; i < e->iw; i++)
	{
		double p[2] = {i, j};
		apply_homography(p, H, p);
		int ip[2] = {p[0], p[1]};
		if (insideP(f, ip[0], ip[1]))
		for (int l = 0; l < 3; l++)
		{
			int idx_win = l + 3 * (f->w * ip[1] + ip[0]);
			float val = getsample_cons(e->img, e->iw, e->ih,
					e->pd, i, j, l);
			f->rgb[idx_win] = 127+val/2;
		}
	}
}

// draw four red segments connecting the control points
static void draw_four_red_segments(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	for (int p = 0; p < 4; p++)
	{
		double x0[2], xf[2];
		map_view_to_window(e, x0, e->c[p]);
		map_view_to_window(e, xf, e->c[(p+1)%4]);
		plot_segment_red(f, x0[0], x0[1], xf[0], xf[1]);
	}
}

// draw disks around the control points
static void draw_four_control_points(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	for (int p = 0; p < 4; p++)
	{
		double P[2];
		map_view_to_window(e, P, e->c[p]);

		// grey circle
		int side = DISK_RADIUS;
		for (int j = -side-1 ; j <= side+1; j++)
		for (int i = -side-1 ; i <= side+1; i++)
		if (hypot(i, j) < side)
		{
			int ii = P[0] + i;
			int jj = P[1] + j;
			if (insideP(f, ii, jj))
				for (int c = 0; c < 3; c++)
					f->rgb[3*(f->w*jj+ii)+c] = 127;
		}

		// central green dot
		int ii = P[0];
		int jj = P[1];
		if (insideP(f, ii, jj))
			f->rgb[3*(f->w*jj+ii)+1]=255;
	}
}

// plot the horizon in green
static void draw_horizon(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	// four control points in projective coordinates
	double a[3] = {e->c[0][0], e->c[0][1], 1};
	double b[3] = {e->c[1][0], e->c[1][1], 1};
	double c[3] = {e->c[2][0], e->c[2][1], 1};
	double d[3] = {e->c[3][0], e->c[3][1], 1};

	// lines though the control points
	double lab[3]; vector_product(lab, a, b);
	double lcd[3]; vector_product(lcd, c, d);
	double lad[3]; vector_product(lad, a, d);
	double lbc[3]; vector_product(lbc, b, c);

	// intersections of opposite sides (vanishing points)
	double p[3]; vector_product(p, lab, lcd);
	double q[3]; vector_product(q, lad, lbc);

	// horizon := line through two vanishing points
	double horizon[3]; vector_product(horizon, p, q);

	// affine coordinates of points
	for (int k = 0; k < 2; k++) {
		p[k] /= p[2];
		q[k] /= q[2];
	}
	p[2] = q[2] = 1;

	// plot the horizon
	double v[2][2];
	map_view_to_window(e, v[0], p);
	map_view_to_window(e, v[1], q);
	if (hypot(hypot(v[0][0], v[0][1]), hypot(v[1][0], v[1][1])) < 1e5)
		plot_segment_green(f, v[0][0], v[0][1], v[1][0], v[1][1]);

}

// Paint the whole scene
// This function is called whenever the window needs to be redisplayed.
static void paint_state(struct FTR *f)
{
	struct viewer_state *e = f->userdata;

	for (int i = 0 ; i < f->w * f->h * 3; i++)
		f->rgb[i] = 255*(i%3); // cyan

	if (e->show_grid_points)
		draw_grid_points(f);
	else
		draw_warped_image(f);
	draw_four_red_segments(f);
	draw_four_control_points(f);
	if (e->show_horizon)
		draw_horizon(f);

	f->changed = 1;
}



// SECTION 8. User-Interface Actions and Events                             {{{1

// action: viewport translation
static void change_view_offset(struct viewer_state *e, double dx, double dy)
{
	e->offset[0] += dx;
	e->offset[1] += dy;
}

// action: viewport zoom
static void change_view_scale(struct viewer_state *e, int x, int y, double fac)
{
	double center[2], X[2] = {x, y};
	map_window_to_view(e, center, X);
	e->scale *= fac;
	for (int p = 0; p < 2; p++)
		e->offset[p] = -center[p]*e->scale + X[p];
	fprintf(stderr, "zoom changed %g\n", e->scale);
}


// test whether (x,y) is inside one of the four control disks
static int hit_point(struct viewer_state *e, double x, double y)
{
	for (int p = 0; p < 4; p++)
	{
		double P[2];
		map_view_to_window(e, P, e->c[p]);
		if (hypot(P[0] - x, P[1] - y) < 2+DISK_RADIUS)
			return p;
	}
	return -1;
}

// key handler
static void event_key(struct FTR *f, int k, int m, int x, int y)
{
	if (k == 'q') {
		ftr_notify_the_desire_to_stop_this_loop(f, 0);
		return;
	}

	struct viewer_state *e = f->userdata;

	if (k == 'c') center_view(f);
	if (k == 'J') change_view_offset(e, 0, -1);
	if (k == 'K') change_view_offset(e, 0, 1);
	if (k == 'H') change_view_offset(e, 1, 0);
	if (k == 'L') change_view_offset(e, -1, 0);
	if (k == 'j') change_view_offset(e, 0, -10);
	if (k == 'k') change_view_offset(e, 0, 10);
	if (k == 'h') change_view_offset(e, 10, 0);
	if (k == 'l') change_view_offset(e, -10, 0);
	if (k == FTR_KEY_DOWN ) change_view_offset(e, 0, -100);
	if (k == FTR_KEY_UP   ) change_view_offset(e, 0, 100);
	if (k == FTR_KEY_RIGHT) change_view_offset(e, -100, 0);
	if (k == FTR_KEY_LEFT)  change_view_offset(e, 100, 0);
	if (k == '+') change_view_scale(e, f->w/2, f->h/2, ZOOM_FACTOR);
	if (k == '-') change_view_scale(e, f->w/2, f->h/2, 1.0/ZOOM_FACTOR);
	if (k == 'p') e->tile_plane = !e->tile_plane;
	if (k == 'w') e->show_horizon = !e->show_horizon;
	if (k >= '0' && k <= '9') e->interpolation_order = k - '0';
	if (k == '.') e->show_grid_points = !e->show_grid_points;
	if (k == 'z') {
		e->dragging_point = false;
		e->dragging_ipoint = false;
	}

	paint_state(f);
}

// resize handler
static void event_resize(struct FTR *f, int k, int m, int x, int y)
{
	paint_state(f);
}

// mouse button handler
static void event_button(struct FTR *f, int k, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// begin dragging a control point in the WINDOW DOMAIN
	if (k == FTR_BUTTON_LEFT)
	{
		int p = hit_point(e, x, y);
		if (p >= 0)
		{
			e->dragging_point = true;
			e->dragged_point = p;
		}
	}

	// end dragging a control point in the WINDOW DOMAIN
	if (e->dragging_point && k == -FTR_BUTTON_LEFT)
	{
		int p = e->dragged_point;
		e->dragging_point = false;
		double X[2] = {x, y};
		map_window_to_view(e, e->c[p], X);
	}

	// begin dragging a control point in the IMAGE DOMAIN
	if (k == FTR_BUTTON_RIGHT)
	{
		int p = hit_point(e, x, y);
		if (p >= 0)
		{
			e->dragging_ipoint = true;
			e->dragged_point = p;
		}
	}

	// end dragging a control point in the IMAGE DOMAIN
	if (e->dragging_ipoint && k == -FTR_BUTTON_RIGHT)
	{
		int p = e->dragged_point;
		e->dragging_ipoint = false;
		double P[2], Q[2] = {x, y};
		map_window_to_image(e, P, Q);
		e->p[p][0] = P[0];
		e->p[p][1] = P[1];
		map_window_to_view(e, e->c[p], Q);
	}

	// zoom in/out
	if (k == FTR_BUTTON_DOWN) change_view_scale(e, x, y, ZOOM_FACTOR);
	if (k == FTR_BUTTON_UP) change_view_scale(e, x, y, 1.0/ZOOM_FACTOR);

	paint_state(f);
}

// mouse motion handler
static void event_motion(struct FTR *f, int b, int m, int x, int y)
{
	struct viewer_state *e = f->userdata;

	// drag WINDOW DOMAIN control point (realtime feedback)
	if (e->dragging_point && m == FTR_BUTTON_LEFT)
	{
		int p = e->dragged_point;
		double X[2] = {x, y};
		map_window_to_view(e, e->c[p], X);
		paint_state(f);
	}

	// drag IMAGE DOMAIN control point (realtime feedback)
	if (e->dragging_ipoint && m == FTR_BUTTON_RIGHT)
	{
		int p = e->dragged_point;
		double P[2], Q[2] = {x, y};
		map_window_to_image(e, P, Q);
		e->p[p][0] = P[0];
		e->p[p][1] = P[1];
		map_window_to_view(e, e->c[p], Q);
		paint_state(f);
	}
}



// SECTION 9. Main Program                                                  {{{1

// library for image input-output
#include "iio.h"

// main function
int main(int argc, char *argv[])
{
	if (argc != 2) {
		fprintf(stderr, "usage:\n\t%s [image.png]\n", *argv);
		return 1;
	}
	char *filename_in = argv[1];

	struct FTR f = ftr_new_window(512,512);
	struct viewer_state e[1];
	f.userdata = e;

	e->img = iio_read_image_float_vec(filename_in, &e->iw, &e->ih, &e->pd);
    float r[(e->iw)*(e->iw)*6+3];
    for(int g = 0;g<(e->iw)*(e->ih)*6;g++){r[g]=0;}
    e->img = mip_map(e->img,e->iw,e->ih,e->pd,r);
    
//	printf("%d %d",e->iw,e->ih);
	center_view(&f);
	paint_state(&f);

	ftr_set_handler(&f, "key", event_key);
	ftr_set_handler(&f, "button", event_button);
	ftr_set_handler(&f, "motion", event_motion);
	ftr_set_handler(&f, "resize", event_resize);

	return ftr_loop_run(&f);
}

// vim:set foldmethod=marker:
