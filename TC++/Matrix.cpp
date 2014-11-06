#include <math.h>

void CrossProduct(float* k,float* a,float* b){
    k[0] = a[1] * b[2] - a[2] * b[1];
    k[1] = a[2] * b[0] - a[0] * b[2];
    k[2] = a[0] * b[1] - a[1] * b[0];
}

float PlaneCrossProduct(float* a,float* b,int k = 0){
    float kk;
	if(k==0){
		kk = a[1]*b[2]-a[2]*b[1];
	}return kk;
}

float DotProduct(float* a,float* b){
	float result;
    result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

    return result;
}

float Modulus(float* a){
	float result;
	result = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

    return result;
}

void Normalize(float *a){
	float Norm = Modulus(a);
	for(int i=0;i<3;i++)	a[i] /= Norm;
}

float NomalizeVector (float* AB,float* A,float* B){
	for (int i=0;i<3;i++){
		AB[i] = B[i] - A[i];
	}
	float norm = Modulus(AB);
	for (int i=0;i<3;i++){
		AB[i] = AB[i]/norm;
	}
	return norm;
}