#include <stdio.h>
#include <math.h>

// Files
extern FILE* log_transform;

// Matrics
extern float rotatinMatrix[9];
extern float translateMatrix[3];

// Matrix.cpp
extern float Normal(float* a);
extern float NomalizeVector (float*,float*,float*);
extern void  CrossProduct(float*,float*,float*);
extern float DotProduct(float* a,float* b);

// Rotate.cpp
extern void  CalcRotate(float*,float*);
extern void  CalcTranslate(float*);
extern void  CalcCoordinate(float*,float*,float*,float*);
extern void  Rodrigues(float cos_theta, float* k);

// Judge.cpp
extern float PlaneCrossProduct(float*,float*,int k = 0);
extern int   JudgePointInQuad(float*,float*,float*,float*);
extern float JudgePointInLine(float*,float*);
extern float JudgePointInTri(float* A,float* B,float* C);

// Interpolate.cpp
extern void  QuadInterFunc(float*,float*,float (*)[3],int k=0);
extern void  TriInterFunc(float *N,float *A,float *B,float *C,int k=0);
extern float OriginPointInterpolateByTri(float*,float*,float*,float*);
extern float Interpolate(float *coor_old,float *N,int k=4);
extern void  CoorInterpolateByTri(float*,float*,float*,float*,float*);