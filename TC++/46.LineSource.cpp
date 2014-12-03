#include "CONSTRUCT.h"
extern float Rb_p, beta_p;

void LineSource(int en, int CrossAxis){
	// 4个积分点方案
	float LMN_Pipe_Center[3], pLMN[3], tmp, PI = 3.1415926f, r = Rb_p/100.0f;
	int num = NodeNum_e[en], m=0, ij[2], k=0;

	float L = 1.0f; // 管长
	float fan = (2*PI) * (r*1.2f) / 4.0f ;	// 1.2倍管周长
	float ww = L * fan;
	float det_r = (r*1.2f)*2.0f/1.5f; // 积分点位置转换到局部坐标下
	GetEleXYZ(8, en);	 // 获取结点坐标

	for(int i=0; i<3; i++)
		if(i!=CrossAxis)	ij[k++] = i;
	
	for(int i=0; i<3; i++)
		LMN_Pipe_Center[i] = 0.5*(LMN_Pipe_e[en][0][i] + LMN_Pipe_e[en][1][i]);

	for(int i=0; i<4; i++){
		fprintf(DEBUG,"%d\n",i);
		float XYZ[3]={0,0,0}, sf_[16], r, dir[3];
		pLMN[ij[0]] = LMN_Pipe_Center[ij[0]] +((i-2)%2)*det_r; // 预定管结点的局部坐标, 以后还要加上由整体坐标获得局部坐标的程序 zuozuoz
		pLMN[ij[1]] = LMN_Pipe_Center[ij[1]] +((i-1)%2)*det_r;
		pLMN[CrossAxis] = LMN_Pipe_Center[CrossAxis];
		fprintf(DEBUG,"=%d\n",i);

		GetXYZ(pLMN, nodeXYZ_e, XYZ);	// 积分点全局坐标
		circular(XYZ, nodeXYZ_e[8], nodeXYZ_e[9], &r, dir);	// 积分点在柱坐标系下的坐标
		Enrichsf(pLMN, r, sf_);

		fprintf(log_bound,"水管积分点 形函数 sf:\n");
		for (int l=0;l<num;l++)	fprintf(log_bound,"%f\t",sf_[l]);
		fprintf(log_bound,"\n水管积分点位置 x=%f, y=%f, z=%f, r=%f, dir=[%f,%f,%f]\n",pLMN[0],pLMN[1],pLMN[2],r,dir[0],dir[1],dir[2]);

		fprintf(log_bound,"\n水管积分点 H矩阵:\n");
		for(int l=0;l<num;l++){
			Q3[en][l] += sf_[l]* ww *diffusivity_m[m][0]/conduction_m[m][0]; // 插值函数 N * 面积微元 / (c * rou); 差一个beta
			for(int r=0;r<num;r++){
				tmp = sf_[r] * sf_[l]* (ww) * beta_p * diffusivity_m[m][0]/conduction_m[m][0];	// 无需利用det转换到全局了，因为本身就是在全局坐标系下讨论
				H[en][l][r] += tmp;
				fprintf(log_bound,"%f\t",tmp);
			}fprintf(log_bound,"\n");
		}
		fprintf(DEBUG,"=%d\n",i);
	}
}