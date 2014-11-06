#include "Data.h"

//������Ԫ�����������������det,sf,dsf,dxyzsfΪ��Ԫ�������	
extern int ijk[3];
extern float nodeXYZ_e[10][3], **det, ***sf, ****dsf, ****dxyzsf;
extern float w2x2x2[8],w3x3x3[27],w4x4x4[64],w5x5x5[125];
extern float i2x2x2_[8][3],i3x3x3_[27][3],i4x4x4_[64][3],i5x5x5_[125][3];

//������Ԫ��������
void GetEleXYZ(int nodenum, int en);
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

//�����κ���log
extern FILE *log_N;