#include "Construct.h"

void HPQ_HEX8(int en, int p, int m, float *w, float **det){
	for(int k=0;k<8;k++){
		Q0[en][k] += rise_m[m]*m_m[m]*sf[en][p][k]*det[en][p]*w[p];//水化热（不含时间项）
		for(int l=0;l<8;l++){
			H[en][k][l] += (diffusivity_m[m][0]*dxyzsf[en][p][0][k]*dxyzsf[en][p][0][l]
				+diffusivity_m[m][1]*dxyzsf[en][p][1][k]*dxyzsf[en][p][1][l]
				+diffusivity_m[m][2]*dxyzsf[en][p][2][k]*dxyzsf[en][p][2][l])
				*det[en][p]*w[p];	// 体积微元需要乘以det，见《有限单元法》 P134
			P[en][k][l]+=sf[en][p][k]*sf[en][p][l]*det[en][p]*w[p];
		}
	}
}

void StiffMatrix_HEX8(int en){
	int m = material_e[en];	// 材料号

	Q0[en] = (float*)calloc(8,sizeof(float));	Alloc2DArray_float(&H[en],8,8);	
	Q3[en] = (float*)calloc(8,sizeof(float));	Alloc2DArray_float(&P[en],8,8);		//PP[i] = Alloc2DArray_float(8,8);

	for(int j=0;j<8;j++){ 
		Q3[en][j]=0;	Q0[en][j]=0;
		for(int k=0;k<8;k++) {	H[en][j][k]=0.0; P[en][j][k]=0.0; }
	}

	for(int j=0;j<(PointNum_e[en]);j++){
		if(plan_e[en]==2) HPQ_HEX8(en, j, m, w2x2x2, det);
		if(plan_e[en]==3) HPQ_HEX8(en, j, m, w3x3x3, det);
		if(plan_e[en]==4) HPQ_HEX8(en, j, m, w4x4x4, det);
		if(plan_e[en]==5) HPQ_HEX8(en, j, m, w5x5x5, det);
	}
}

void HPQ_HEX16(int en, int p, int m, float *w, float *det){
	for (int j=0;j<16;j++){	// 结点数
		for (int k=0;k<16;k++){
			H[en][j][k] +=	(diffusivity_m[m][0]*dxyzsf[en][p][0][j]*dxyzsf[en][p][0][k]	// 偏x
							+diffusivity_m[m][1]*dxyzsf[en][p][1][j]*dxyzsf[en][p][1][k]	// 偏y
							+diffusivity_m[m][2]*dxyzsf[en][p][2][j]*dxyzsf[en][p][2][k])	// 偏z
							*det[p]*w[p];	// 面积微元×积分权系数， 见《有限单元法》P135
			PP[en][j][k] += sf[en][p][j]*sf[en][p][k]*det[p]*w[p];
		}
	}
}

void StiffMatrix_HEX16(int en){
	int m = material_e[en];	// 材料号
	Q0[en] = (float*)calloc(16,sizeof(float));	Alloc2DArray_float(&H[en],16,16);
	Q3[en] = (float*)calloc(16,sizeof(float));	Alloc2DArray_float(&P[en],16,16);	Alloc2DArray_float(&PP[en],16,16);

	for(int j=0;j<16;j++){
		Q3[en][j]=0;		Q0[en][j]=0;
		for(int k=0;k<8;k++){
			H[en][j][k]=0.0;		P[en][j][k]=0.0;		PP[en][j][k]=0.0;
		}
	}
	for (int j=0;j<PointNum_e[en];j++){
		for(int k=0;k<3;k++){
			for(int l=0;l<8;l++){
				fprintf(DEBUG,"dxyzdsf[%d][%d][%d]=%f \n",j,k,l,dxyzsf[en][j][k][l]);
			}
		}
	}
	for(int i=0;i<PointNum_e[en];i++){
		if(plan_e[en]==2) HPQ_HEX16(en, i, m, w[en], det[en]);
		if(plan_e[en]==3) HPQ_HEX16(en, i, m, w[en], det[en]);
		if(plan_e[en]==4) HPQ_HEX16(en, i, m, w[en], det[en]);
		if(plan_e[en]==5) HPQ_HEX16(en, i, m, w[en], det[en]);
	}
}