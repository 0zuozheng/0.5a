#include "ELEMENT.h"
extern float Rb_p;
extern float Modulus(float*);
extern void Normalize(float*);

void circular(float XYZ[3],float P1[3],float P2[3],float *r,float d[3]){ // 一定要在同一坐标系下，即水管结点和积分点都在全局坐标系
	float PI = 3.1415926;	float XYZ_0[3]; float sum=0.0;
	float gx = XYZ[0]-P1[0],  gy = XYZ[1]-P1[1], gz = XYZ[2]-P1[2];
	float lx = P2[0]-P1[0],   ly = P2[1]-P1[1],  lz = P2[2]-P1[2];
	float A = (lx*gx+ly*gy+lz*gz)/(lx*lx+ly*ly+lz*lz);
	XYZ_0[0] = A*lx + P1[0];
	XYZ_0[1] = A*ly + P1[1];
	XYZ_0[2] = A*lz + P1[2];
	for(int i=0;i<3;i++)
		d[i] = XYZ[i] - XYZ_0[i];
	*r = Modulus(d);
	Normalize(d);
	*r /= Rb_p;
}

void Enrichsf(float xyz[3],float r,float *sf_){ // 丰富形函数
	GetSf(xyz[0],xyz[1],xyz[2],sf_);	// 前8项

	float b = 6.4f;	// ZUOZUOZ
	float phi = - exp(-b*r); // f(r)=-exp(-b*r), b=6.4
	
	for (int i=0;i<8;i++)	sf_[i+8] = phi*sf_[i];	// 扩充形函数
}

void Enrichdxyzdsf(float r,float *dir,float *sf_,float **dsf_,float *det_,float **dxydsf_){
	float **dxydsf_t, dphi[3];
	Alloc2DArray_float(&dxydsf_t,3,8);
	GetDetDxyzsf(det_,dsf_,dxydsf_t);

	float b=6.4;
	float phi = - exp(-b*r);

	extern float Rb_p;

	dphi[0] = b*exp(-b*r)*dir[0]/Rb_p;
	dphi[1] = b*exp(-b*r)*dir[1]/Rb_p;
	dphi[2] = b*exp(-b*r)*dir[2]/Rb_p;

	for (int i=0;i<3;i++){
		for(int j=0;j<8;j++){
			dxydsf_[i][j] = dxydsf_t[i][j];
			dxydsf_[i][j+8] = phi*dxydsf_[i][j] + dphi[i]*sf_[j];
		}
	}
}