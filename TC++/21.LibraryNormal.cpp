/* * * * * * * * * * * * * * * * * *
*
*		常规单元库
*
* * * * * * * * * * * * * * * * * */
#include "ELEMENT.h"
extern void N_Hammer(int en, int order, int GaussNum, float Qi[][3], int CrossFace);

void write_element(int en, int in, int classnum){
	fprintf(log_N,"\n\tGauss Point [%d]\n",in);
	float xx=0,yy=0,zz=0;
	for (int k=0;k<8;k++){
		xx += sf[en][in][k]*nodeXYZ_e[k][0];	yy += sf[en][in][k]*nodeXYZ_e[k][1];	zz += sf[en][in][k]*nodeXYZ_e[k][2];
	}fprintf(log_N,"\t\tLocation：x = %f, y = %f, z = %f\n",xx,yy,zz);

	fprintf(log_N,"\t\tShapeFunc sf [1~%d]:\n",classnum);
	for(int k=0;k<classnum;k++)	fprintf(log_N,"%f\t",sf[en][in][k]);
	fprintf(log_N,"\n\t\tDerivative: dsf\n");
	for (int k=0;k<3;k++){
		for (int l=0;l<classnum;l++)	fprintf(log_N,"%f\t",dsf[en][in][k][l]);
		fprintf(log_N,"\n");
	}fprintf(log_N,"\t\tDerivative: dxydsf\n");
	for (int k=0;k<3;k++){
		for (int l=0;l<classnum;l++)	fprintf(log_N,"%f\t",dxyzsf[en][in][k][l]);
		fprintf(log_N,"\n");
	}
}

void N_Gauss(int en, int order, int GaussNum, float Qi[][3]){
	fprintf(log_N,"\nElementID: [%d] Class: %d Gauss Type: %d\n{",en,class_e[en],plan_e[en]);
	GetEleXYZ(8,en);	// 获取结点坐标

	det[en] = (float *)calloc(GaussNum,sizeof(float));		Alloc2DArray_float(&sf[en],GaussNum,8);
	Alloc3DArray_float(&dsf[en],GaussNum,3,8);				Alloc3DArray_float(&dxyzsf[en],GaussNum,3,8);

	for(int j=0;j<GaussNum;j++){
		GetSf(Qi[j][0], Qi[j][1], Qi[j][2], sf[en][j]);
		GetDsf(Qi[j][0], Qi[j][1], Qi[j][2], dsf[en][j]);
		GetDetDxyzsf(&det[en][j], dsf[en][j], dxyzsf[en][j]);
		write_element(en,j,8);
	}fprintf(log_N,"}\n");	for(int i=0;i<90;i++)	fprintf(log_N,"-");
}

void ShapeFunc_normal(int en){
	if (plan_e[en]==1008)	N_Gauss(en,0, 8,i2x2x2_);
	if (plan_e[en]==1027)	N_Gauss(en,0,27,i3x3x3_);
	if (plan_e[en]==1064)	N_Gauss(en,0,64,i4x4x4_);
	if (plan_e[en]==1125)	N_Gauss(en,0,125,i5x5x5_);
}

void FreeShape(){
	int GaussNum=8;
	for(int i=0;i<elementnum_c;i++){
		GaussNum = plan_e[i]%1000;
		free(det[i]);
		FREE2D_float(sf[i],GaussNum);
		FREE3D_float(dsf[i],GaussNum,3);
		FREE3D_float(dxyzsf[i],GaussNum,3);
	}
	free(det);	free(sf);
	free(dsf);	free(dxyzsf);
}