/* * * * * * * * * * * * * * * * * *
*
*		常规单元库
*
* * * * * * * * * * * * * * * * * */
#include "ELEMENT.h"

void write_element(int en, int in, int classnum){
	fprintf(log_N,"\n积分点 [%d]\n",in);
	float xx=0,yy=0,zz=0;
	for (int k=0;k<8;k++){
		xx += sf[en][in][k]*nodexyz_e[0][k];	yy += sf[en][in][k]*nodexyz_e[1][k];	zz += sf[en][in][k]*nodexyz_e[2][k];
	}fprintf(log_N," 积分点位置：\n x = %f y = %f z = %f\n",xx,yy,zz);

	fprintf(log_N," 形函数 sf [1~%d]:\n",classnum);
	for(int k=0;k<classnum;k++)	fprintf(log_N,"%f\t",sf[en][in][k]);
	fprintf(log_N,"\n 偏导数 dsf\n");
	for (int k=0;k<3;k++){
		for (int l=0;l<classnum;l++)	fprintf(log_N,"%f\t",dsf[en][in][k][l]);
		fprintf(log_N,"\n");
	}fprintf(log_N," 偏导数 dxydsf\n");
	for (int k=0;k<3;k++){
		for (int l=0;l<classnum;l++)	fprintf(log_N,"%f\t",dxyzsf[en][in][k][l]);
		fprintf(log_N,"\n");
	}
}

void N_1(int en, int order, int GaussNum, float Qi[][3]){
	fprintf(log_N,"\n单元 - ElementID: [%d] 单元类型 - Class: [%d]",en,class_e[en]);
	GetNodeXYZ(8,en);

	det[en] = (float *)calloc(GaussNum,sizeof(float));		Alloc2DArray_float(&sf[en],GaussNum,8);
	Alloc3DArray_float(&dsf[en],GaussNum,3,8);				Alloc3DArray_float(&dxyzsf[en],GaussNum,3,8);

	for(int j=0;j<GaussNum;j++){
		GetSfDsf(Qi[j][0], Qi[j][1], Qi[j][2], sf[en][j], dsf[en][j]);
		GetDetDxyzsf(&det[en][j], dsf[en][j], dxyzsf[en][j]);

		write_element(en,j,8);
	}
}

void ShapeFunc_normal(int en){
	if (plan_e[en]==1008)	N_1(en,0, 8,i2x2x2_);
	if (plan_e[en]==1027)	N_1(en,0,27,i3x3x3_);
	if (plan_e[en]==1064)	N_1(en,0,64,i4x4x4_);
	if (plan_e[en]==1125)	N_1(en,0,125,i5x5x5_);
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
	free(det);
	free(sf);
	free(dsf);
	free(dxyzsf);
}