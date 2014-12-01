/* * * * * * * * * * * * * * * * * *
*
*		单元参数
*
* * * * * * * * * * * * * * * * * */
#include <stdio.h>
extern FILE *DEBUG;

int   ijk[3] = {1,2,0};
float w2x2x2[8],w3x3x3[27],w4x4x4[64],w5x5x5[125];
float i2x2x2_[8][3],i3x3x3_[27][3],i4x4x4_[64][3],i5x5x5_[125][3];
float tw2xw2[6],   tw3xw3[12],   tw4xw4[24],   tw5xw5[35],    tw6xw6[72],   tw7xw7[91];
float ti2xi2[6][4],ti3xi3[12][4],ti4xi4[24][4],ti5xi5[35][4], ti6xi6[72][4],ti7xw7[91][4];

float w2[2] = {     1.0f,    1.0f};	// 积分权系数(2个高斯积分点)
float i2[2] = {-0.57735f,0.57735f};
float w3[3] = { 0.55555f, 0.88888f, 0.55555f};
float i3[3] = {-0.77460f, 0.0f,     0.77460f};
float w4[4] = { 0.34785f, 0.65214f, 0.65214f, 0.34785f};
float i4[4] = {-0.86113f,-0.33998f, 0.33998f, 0.86113f};
float w5[5] = {0.236927f, 0.47863f, 0.56889f, 0.47863f, 0.236927f};
float i5[5] = {-0.90618f,-0.53847f, 0.0f,     0.53847f,  0.90618f};

int PointPerTri[21] = {0, 1, 3, 4, 6, 7, 12, 13, 16, 19, 25, 27, 33, 37, 42, 48, 52, 61, 70, 73, 79};	// 各阶精度对应的积分点数 直到20阶次
float tw2[3] = {0.33333f,0.33333f,0.33333f};
float ti2[3][3] = {{ 0.5f, 0.5f, 0.0f},
					{ 0.0f, 0.5f, 0.5f},
					{ 0.5f, 0.0f, 0.5f}};
float tw3[4] = {-0.5625f,0.5208f,0.5208f,0.5208f};
float ti3[4][3] = {{ 0.33333f, 0.33333f, 0.33333f},
					{ 0.73333f, 0.13333f, 0.13333f},
					{ 0.13333f, 0.73333f, 0.13333f},
					{ 0.13333f, 0.13333f, 0.73333f}};
float tw4[6]={0.22338f,0.22338f,0.22338f,0.10995f,0.10995f,0.10995f};
float ti4[6][3] = {{ 0.10810f,  0.44595f,  0.44595f},
					{ 0.44595f,  0.10810f,  0.44595f},
					{ 0.44595f,  0.44595f,  0.10810f},
					{ 0.81684f,  0.09158f,  0.09158f},
					{ 0.09158f,  0.81684f,  0.09158f},
					{ 0.09158f,  0.09158f,  0.81684f}};
float tw5[7]={0.225f,0.1324f,0.1324f,0.1324f,0.12594f,0.12594f,0.12594f};
float ti5[7][3] = {{ 0.33333f,  0.33333f,  0.33333f},
					{ 0.05971f, 0.47014f,  0.47014f},
					{ 0.47014f, 0.05971f,  0.47014f},
					{ 0.47014f, 0.47014f,  0.05971f},
					{ 0.797427f, 0.10129f,  0.10129f},
					{ 0.10129f,  0.797427f, 0.10129f},
					{ 0.10129f,  0.10129f,  0.797427f}};
float tw7[13]={-0.14957f,0.17561f,0.17561f,0.17561f,0.05335f,0.05335f,0.05335f,0.07711f,0.07711f,0.07711f,0.07711f,0.07711f,0.07711f};
float ti7[13][3]={{ 0.33333f,  0.33333f,  0.33333f},
					{ 0.47931f, 0.260345f, 0.260345f},
					{0.260345f,  0.47931f, 0.260345f},
					{0.260345f, 0.260345f,  0.47931f},
					{ 0.86974f,  0.06513f,  0.06513f},
					{ 0.06513f,  0.86974f,  0.06513f},
					{ 0.06513f,  0.06513f,  0.86974f},
					{ 0.63844f,  0.31286f,  0.04869f},
					{ 0.63844f,  0.04869f,  0.31286f},
					{ 0.31286f,  0.63844f,  0.04869f},
					{ 0.04869f,  0.63844f,  0.31286f},
					{ 0.31286f,  0.04869f,  0.63844f},
					{ 0.04869f,  0.31286f,  0.63844f}};

extern int Axis_Of_Face[65];

void IP_HEX_Gauss(int Num,float* ii,float* w,float ixxx[][3],float* wxxx){
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

void IP_Prism_Hammer(int order,float *i_line,float *w_line,float i_tri[][3],float *w_tri,float ixxx[][4],float *wxxx){
	int TriNum = PointPerTri[order];
	for (int i=0;i<order;i++){
		for (int j=0;j<TriNum;j++){
			wxxx[TriNum * i + j]    = w_line[i] * w_tri[j];
			ixxx[TriNum * i + j][0] = i_tri[j][0];
			ixxx[TriNum * i + j][1] = i_tri[j][1];
			ixxx[TriNum * i + j][2] = i_tri[j][2];
			ixxx[TriNum * i + j][3] = i_line[i];
		}
	}
}

void IntegratePlan(){ // 对各种积分方案的常数预先处理
	IP_HEX_Gauss(2,i2,w2,i2x2x2_,w2x2x2);
	IP_HEX_Gauss(3,i3,w3,i3x3x3_,w3x3x3);
	IP_HEX_Gauss(4,i4,w4,i4x4x4_,w4x4x4);
	IP_HEX_Gauss(5,i5,w5,i5x5x5_,w5x5x5);

	IP_Prism_Hammer(2,i2,w2,ti2,tw2,ti2xi2,tw2xw2);
	IP_Prism_Hammer(3,i3,w3,ti3,tw3,ti3xi3,tw3xw3);
	IP_Prism_Hammer(3,i4,w4,ti4,tw4,ti4xi4,tw4xw4);
	IP_Prism_Hammer(5,i5,w5,ti5,tw5,ti5xi5,tw5xw5);
	for(int i=0;i<35;i++)
		fprintf(DEBUG,"w = %f, i = [%f %f %f %f]\n",tw5xw5[i],ti5xi5[i][0],ti5xi5[i][1],ti5xi5[i][2],ti5xi5[i][3]);
}