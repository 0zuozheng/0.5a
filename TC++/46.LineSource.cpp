#include "CONSTRUCT.h"
extern float Rb_p, beta_p;

void LineSource(int en, int CrossAxis){
	// 4�����ֵ㷽��
	float LMN_Pipe_Center[3], pLMN[3], tmp, PI = 3.1415926f, r = Rb_p/100.0f;
	int num = NodeNum_e[en], m=0, ij[2], k=0;

	float L = 1.0f; // �ܳ�
	float fan = (2*PI) * (r*1.2f) / 4.0f ;	// 1.2�����ܳ�
	float ww = L * fan;
	float det_r = (r*1.2f)*2.0f/1.5f; // ���ֵ�λ��ת�����ֲ�������
	GetEleXYZ(8, en);	 // ��ȡ�������

	for(int i=0; i<3; i++)
		if(i!=CrossAxis)	ij[k++] = i;
	
	for(int i=0; i<3; i++)
		LMN_Pipe_Center[i] = 0.5*(LMN_Pipe_e[en][0][i] + LMN_Pipe_e[en][1][i]);

	for(int i=0; i<4; i++){
		fprintf(DEBUG,"%d\n",i);
		float XYZ[3]={0,0,0}, sf_[16], r, dir[3];
		pLMN[ij[0]] = LMN_Pipe_Center[ij[0]] +((i-2)%2)*det_r; // Ԥ���ܽ��ľֲ�����, �Ժ�Ҫ���������������þֲ�����ĳ��� zuozuoz
		pLMN[ij[1]] = LMN_Pipe_Center[ij[1]] +((i-1)%2)*det_r;
		pLMN[CrossAxis] = LMN_Pipe_Center[CrossAxis];
		fprintf(DEBUG,"=%d\n",i);

		GetXYZ(pLMN, nodeXYZ_e, XYZ);	// ���ֵ�ȫ������
		circular(XYZ, nodeXYZ_e[8], nodeXYZ_e[9], &r, dir);	// ���ֵ���������ϵ�µ�����
		Enrichsf(pLMN, r, sf_);

		fprintf(log_bound,"ˮ�ܻ��ֵ� �κ��� sf:\n");
		for (int l=0;l<num;l++)	fprintf(log_bound,"%f\t",sf_[l]);
		fprintf(log_bound,"\nˮ�ܻ��ֵ�λ�� x=%f, y=%f, z=%f, r=%f, dir=[%f,%f,%f]\n",pLMN[0],pLMN[1],pLMN[2],r,dir[0],dir[1],dir[2]);

		fprintf(log_bound,"\nˮ�ܻ��ֵ� H����:\n");
		for(int l=0;l<num;l++){
			Q3[en][l] += sf_[l]* ww *diffusivity_m[m][0]/conduction_m[m][0]; // ��ֵ���� N * ���΢Ԫ / (c * rou); ��һ��beta
			for(int r=0;r<num;r++){
				tmp = sf_[r] * sf_[l]* (ww) * beta_p * diffusivity_m[m][0]/conduction_m[m][0];	// ��������detת����ȫ���ˣ���Ϊ���������ȫ������ϵ������
				H[en][l][r] += tmp;
				fprintf(log_bound,"%f\t",tmp);
			}fprintf(log_bound,"\n");
		}
		fprintf(DEBUG,"=%d\n",i);
	}
}