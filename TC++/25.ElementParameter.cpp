/* * * * * * * * * * * * * * * * * *
*
*		单元参数
*
* * * * * * * * * * * * * * * * * */
#include <stdio.h>

int ijk[3]={1,2,0};
float w2x2x2[8],w3x3x3[27],w4x4x4[64],w5x5x5[125];
float i2x2x2_[8][3],i3x3x3_[27][3],i4x4x4_[64][3],i5x5x5_[125][3];

float w2[2] = {     1.0f,    1.0f};	// 积分权系数(2个高斯积分点)
float i2[2] = {-0.57735f,0.57735f};
float w3[3] = { 0.55555f, 0.88888f, 0.55555f};
float i3[3] = {-0.77460f, 0.0f,     0.77460f};
float w4[4] = { 0.34785f, 0.65214f, 0.65214f, 0.34785f};
float i4[4] = {-0.86113f,-0.33998f, 0.33998f, 0.86113f};
float w5[5] = {0.236927f, 0.47863f, 0.56889f, 0.47863f, 0.236927f};
float i5[5] = {-0.90618f,-0.53847f, 0.0f,     0.53847f,  0.90618f};

void IntPoint_Gauss(int Num,float* ii,float* w,float ixxx[][3],float* wxxx){
	for (int i=0;i<Num;i++){
		for (int j=0;j<Num;j++){
			for (int k=0;k<Num;k++){
				wxxx[Num*Num*i+Num*j+k] = w[i] * w[j] * w[k];
				ixxx[Num*Num*i+Num*j+k][0] = ii[i];
				ixxx[Num*Num*i+Num*j+k][1] = ii[j];
				ixxx[Num*Num*i+Num*j+k][2] = ii[k];
			}
		}
	}
}

void IntegratePlan(){
	IntPoint_Gauss(2,i2,w2,i2x2x2_,w2x2x2);
	IntPoint_Gauss(3,i3,w3,i3x3x3_,w3x3x3);
	IntPoint_Gauss(4,i4,w4,i4x4x4_,w4x4x4);
	IntPoint_Gauss(5,i5,w5,i5x5x5_,w5x5x5);
}