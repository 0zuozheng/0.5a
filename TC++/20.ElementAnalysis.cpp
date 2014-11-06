/* * * * * * * * * * * * * * * * * *
*
*		单元推导
*
* * * * * * * * * * * * * * * * * */

#include "ELEMENT.h"

float nodeXYZ_e[10][3], jacobi[3][3], divjacobi[3][3], **det, ***sf, ****dsf, ****dxyzsf;

void GetEleXYZ(int nodenum,int en){ //en -> element_now;	该子程序只计算当前单元的Node信息，每次调用都重置一遍
	for(int i=0;i<nodenum;i++){
		int k = node_e[en][i];	//获取结点号
		for(int j=0;j<3;j++)	nodeXYZ_e[i][j]=xyz_n[k][j];	//获取结点坐标，赋予nodexyz_e
	}
}

void GetSf(float x, float y, float z, float* sf_){ //x,y,z -> coordinate; sf_ -> 形函数; dsf -> 形函数的偏导数(局部坐标系),即 (偏N/偏epsilon).
	for(int i=0;i<8;i++){
		sf_[i]=0;
	}	// Serendipity积分, 《有限单元法》 P114
	sf_[0]=(1-x)*(1-y)*(1-z)/8;	 sf_[1]=(1+x)*(1-y)*(1-z)/8;    sf_[2]=(1+x)*(1+y)*(1-z)/8;	 sf_[3]=(1-x)*(1+y)*(1-z)/8;
	sf_[4]=(1-x)*(1-y)*(1+z)/8;	 sf_[5]=(1+x)*(1-y)*(1+z)/8;    sf_[6]=(1+x)*(1+y)*(1+z)/8;	 sf_[7]=(1-x)*(1+y)*(1+z)/8;
}	//以上子程序中均在局部坐标系下讨论

void GetDsf(float x, float y, float z, float **dsf_){ //x,y,z -> coordinate; sf_ -> 形函数; dsf -> 形函数的偏导数(局部坐标系),即 (偏N/偏epsilon).
	for(int i=0;i<8;i++){
		for(int j=0;j<3;j++)  dsf_[j][i]=0;	//矩阵置零
	}	// Serendipity积分, 《有限单元法》 P114
	dsf_[0][0]=-(1-y)*(1-z)/8;	 dsf_[0][1]=(1-y)*(1-z)/8;     dsf_[0][2]=(1+y)*(1-z)/8;	  dsf_[0][3]=-(1+y)*(1-z)/8;
	dsf_[0][4]=-(1-y)*(1+z)/8;	 dsf_[0][5]=(1-y)*(1+z)/8;     dsf_[0][6]=(1+y)*(1+z)/8;	  dsf_[0][7]=-(1+y)*(1+z)/8;

	dsf_[1][0]=-(1-x)*(1-z)/8;	 dsf_[1][1]=-(1+x)*(1-z)/8;    dsf_[1][2]=(1+x)*(1-z)/8;	  dsf_[1][3]=(1-x)*(1-z)/8;
	dsf_[1][4]=-(1-x)*(1+z)/8;	 dsf_[1][5]=-(1+x)*(1+z)/8;    dsf_[1][6]=(1+x)*(1+z)/8;	  dsf_[1][7]=(1-x)*(1+z)/8;

	dsf_[2][0]=-(1-x)*(1-y)/8;	 dsf_[2][1]=-(1+x)*(1-y)/8;	  dsf_[2][2]=-(1+x)*(1+y)/8;	  dsf_[2][3]=-(1-x)*(1+y)/8;
	dsf_[2][4]=(1-x)*(1-y)/8;	 dsf_[2][5]=(1+x)*(1-y)/8;	  dsf_[2][6]=(1+x)*(1+y)/8;	      dsf_[2][7]=(1-x)*(1+y)/8;
}	//以上子程序中均在局部坐标系下讨论

void GetDetDxyzsf(float *det_, float **dsf_, float **dxyzsf_){ //det ->	Jacobi阵的行列式	dxyzsf -> 形函数偏导(整体坐标系)
	int i,j,k;	float tmp;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++)  { jacobi[j][i]=0; divjacobi[j][i]=0; }// 重置Jacobi阵
		for(j=0;j<8;j++)  dxyzsf_[i][j]=0;						// 重置形函数偏导阵(整体坐标系)
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tmp=0.0;
			for(k=0;k<8;k++)  tmp+=dsf_[i][k]*nodeXYZ_e[k][j];	// 这个是每个单元都不同的
			jacobi[i][j]=tmp;	//生成Jacobi阵
		}
	}
	for(i=0;i<3;i++){	//计算伴随矩阵
		j = ijk[i];   k = ijk[j];
		divjacobi[i][i] = jacobi[j][j]*jacobi[k][k]-jacobi[k][j]*jacobi[j][k];
		divjacobi[i][j] = jacobi[k][j]*jacobi[i][k]-jacobi[i][j]*jacobi[k][k];
		divjacobi[j][i] = jacobi[j][k]*jacobi[k][i]-jacobi[j][i]*jacobi[k][k];
	}
	*det_ = jacobi[0][0]*divjacobi[0][0]+jacobi[0][1]*divjacobi[1][0]+jacobi[0][2]*divjacobi[2][0];	//计算行列式
	

	for(i=0;i<3;i++){
		for(j=0;j<8;j++){
			tmp = 0;
			for(k=0;k<3;k++)  tmp += divjacobi[i][k]*dsf_[k][j];
			dxyzsf_[i][j] = tmp/(*det_);	// 逆矩阵 = 伴随矩阵 / 行列式
		}
	}
}

void ElementAnalysis(){	// Called by Main.cpp
	det = (float **)calloc(elementnum_c,sizeof(float *));
	sf  = (float ***)calloc(elementnum_c,sizeof(float **));
	dsf  = (float ****)calloc(elementnum_c,sizeof(float ***));
	dxyzsf  = (float ****)calloc(elementnum_c,sizeof(float ***));
	
	IntegratePlan();	// 对各种积分方案的常数预先处理, 包括i2x2x2, w2x2x2之类的

	for(int i=0;i<elementnum_c;i++)	{
		if(pipe_e[i])	ShapeFunc_enrich(i);
		else			ShapeFunc_normal(i);
	}
	fclose(log_N);
}