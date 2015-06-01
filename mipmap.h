//Mip-Map section by Hugo

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

//pour plus de modularité, on déclare ce filtre qui est la moyenne
//float filter_mipmap(int i,int j,int w,int pd,int d,int l,float *r){
//	return (r[coord_mipmap(2*i,2*j,w,pd,d-1,l)]+r[coord_mipmap(2*i+1,2*j,w,pd,d-1,l)]+r[coord_mipmap(2*i,2*j+1,w,pd,d-1,l)]+r[coord_mipmap(2*i+1,2*j+1,w,pd,d-1,l)])/4;
//}

#define TAPS 7 //racine du nombre de coefficient interpolé, doit être impaire

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


#define SIG 0.6 //paramètre du filtre gaussien

// on admet pour l'instant que l'image et de taille une puissance de 2, cette fonction construit le mip_map d'une image
float *mip_map(float *img,int w,int h, int pd, float *r){
int i,j,d,l,ll;
double gauss[TAPS*TAPS];
for (i=0;i<TAPS;i++){
	for(j=0;j<TAPS;j++){
		gauss[i+TAPS*j]=exp(-(pow(j-(TAPS-1)/2,2)+pow(i-(TAPS-1)/2,2))/(2*pow(SIG,2)))/(2*M_PI*pow(SIG,2));   //on a prit sigma = 0,6 (ptet le prendre variable)
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


void precal_D(double H[3][3],double *D){
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

#define D_BIAS 0.
// elle prend des int est retourne un double
double cal_D(int x,int y,double *D,double coo[2]){
	double p[2]={x,y};
	double a,b,c;
	c = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	//On declare les dérivé partielles	
	double dudx,dudy,dvdx,dvdy;
	dudx = (D[1]*p[1]+D[2])/c;
	dudy = (D[5]*p[0]+D[6])/c;
	dvdx = (D[3]*p[1]+D[4])/c;
	dvdy = (D[7]*p[0]+D[8])/c;
	
//On decale l'origine du rectangle pour que le parallélogramme soit à l'intérieur	
	if(dudx<0){coo[0] += dudx;}
	if(dudy<0){coo[0] += dudy;}
	if(dvdx<0){coo[1] += dvdx;}
	if(dvdy<0){coo[1] += dvdy;}
		
	a = sqrt (pow(dudx,2) + pow(dvdx,2));
	b = sqrt (pow(dudy,2) + pow(dvdy,2));
	if(a>b){return a+D_BIAS;}{return b+D_BIAS;}
}


//int coord_mipmap (int u,int v,int w, int pd,int d,int l)
// for bilinear interpolation during mipmap
static float bilinear_mipmap(float *x, int w, int h,
		float p, float q, int l, extrapolator_t pix, int d)
{
	float pp = p - 1/2*(pow(2,d)-1)/pow(2,d);
	float qq = q - 1/2*(pow(2,d)-1)/pow(2,d);
	int ip = floor(pp);
	int iq = floor(qq);
	float a = x[coord_mipmap(ip,iq,w,3,d,l)];
	float b = x[coord_mipmap(ip+1,iq,w,3,d,l)];
	float c = x[coord_mipmap(ip,iq+1,w,3,d,l)];
	float dd = x[coord_mipmap(ip+1,iq+1,w,3,d,l)];
	return evaluate_bilinear_cell(a, b, c, dd, p-ip, q-iq);
}



//mipmapping interpolator, img is a mipmapping
static float mip_mapping_interpolator_at(float *img, int w, int h, int pd,
		float x, float y, int l, extrapolator_t p,double D){
	int d ;
	float a,b;
	for(d=0;(pow(2,d))<=D && (pow(2,d))<=w;d++){;}
	if (d==1||d==0){return bilinear_mipmap(img,w,h,x,y,l,p,1);}
	if ((pow(2,d))>=w){return img[6*w*(w-1)+l];} 
	a = bilinear_mipmap(img,w,h,x/pow(2,d-1),y/pow(2,d-1),l,p,d);
	b = bilinear_mipmap(img,w,h,x/pow(2,d-2),y/pow(2,d-2),l,p,d-1);
	return ((D-pow(2,d-1))/(pow(2,d-1)))*a + ((pow(2,d)-D)/pow(2,d-1))*b;
}


// on a l'impression que img[i][j][l] avec l la couleur, pd le nombre de couleur = img[(i+j*w)*pd + l]
// et dans le mip-map ?
// r[u][v][d] = h*(w+..+w/(2^(d-1)))*pd + (u+v*w/d)*pd+l


// on note D = (1,2,3,4,5,6,7,8,9,10,11) avec 
// D = max(sqrt( (1y+2)^2 + (3y+4)^2 )),sqrt( (5x+6)^2 + (7x+8)^2 ))/(9x + 10y +11)^2

