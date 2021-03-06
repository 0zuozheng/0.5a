#include "Data.h"

//声明单元分析计算变量，其中det,sf,dsf,dxyzsf为单元分析结果	
extern int ijk[3];
extern float nodeXYZ_e[10][3], **w, **det, ***sf, ****dsf, ****dxyzsf;
extern float w2x2x2[8],w3x3x3[27],w4x4x4[64],w5x5x5[125];
extern float i2x2x2_[8][3],i3x3x3_[27][3],i4x4x4_[64][3],i5x5x5_[125][3];
extern float tw2xw2[6],   tw3xw3[12],   tw4xw4[24],   tw5xw5[35],    tw6xw6[72],   tw7xw7[91];
extern float ti2xi2[6][4],ti3xi3[12][4],ti4xi4[24][4],ti5xi5[35][4], ti6xi6[72][4],ti7xw7[91][4];

//声明单元分析计算
void GetEleXYZ(int nodenum, int en);
void GetXYZ(float coor[3],float (*XYZ_e)[3],float XYZ[3]);
void GetSf(float x, float y, float z, float* sf_);
void GetDsf(float x, float y, float z, float** dsf_);
void GetDetDxyzsf(float* det_, float** dsf_, float** dxyzsf_);
void ElementAnalysis();

//21.LibraryNormal.cpp
extern void ShapeFunc_normal(int);
extern void ShapeFunc_enrich(int);

//25.ElementParameter.cpp
extern void IntegratePlan();
extern void AxisPlan();
extern void PrismOrder(int order,float (*i_tri)[3],float *i_line);

//声明形函数log
extern FILE *log_N;