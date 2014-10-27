#include "Construct.h"

void HPQ_HEX8(int i, int j, int m, float *w, float **det){
	for(int k=0;k<8;k++){
		Q0[i][k]+=rise_m[m]*m_m[m]*sf[i][j][k]*det[i][j]*w[j];//水化热（不含时间项）
		for(int l=0;l<8;l++){
			H[i][k][l]+=(diffusivity_m[m][0]*dxyzsf[i][j][0][k]*dxyzsf[i][j][0][l]
				+diffusivity_m[m][1]*dxyzsf[i][j][1][k]*dxyzsf[i][j][1][l]
				+diffusivity_m[m][2]*dxyzsf[i][j][2][k]*dxyzsf[i][j][2][l])
				*det[i][j]*w[j];	// 体积微元需要乘以det，见《有限单元法》 P134
			P[i][k][l]+=sf[i][j][k]*sf[i][j][l]*det[i][j]*w[j];
		}
	}
}

void StiffMatrix_HEX8(int i){
	int m=material_e[i];	// 材料号

	Q0[i] = (float*)calloc(8,sizeof(float));	Alloc2DArray_float(&H[i],8,8);	
	Q3[i] = (float*)calloc(8,sizeof(float));	Alloc2DArray_float(&P[i],8,8);		//PP[i] = Alloc2DArray_float(8,8);

	for(int j=0;j<8;j++){ 
		Q3[i][j]=0;	Q0[i][j]=0;
		for(int k=0;k<8;k++) {	H[i][j][k]=0.0; P[i][j][k]=0.0; }
	}

	for(int j=0;j<(plan_e[i]%1000);j++){
		if(plan_e[i]==1008) HPQ_HEX8(i, j, m, w2x2x2, det);
		if(plan_e[i]==1027) HPQ_HEX8(i, j, m, w3x3x3, det);
		if(plan_e[i]==1064) HPQ_HEX8(i, j, m, w4x4x4, det);
		if(plan_e[i]==1125) HPQ_HEX8(i, j, m, w5x5x5, det);
	}
}