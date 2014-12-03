#include "Construct.h"
float beta_p = 50663.0f;		// lmd = 1, beta = lmd/(1.2*r0*ln(r0/rp)), 185/(1.2*0.84/100*log(1.2))
float Rb_p   = 0.84f;			// 对应1.5m×1.5m

void pipe_alongT(int PipeId, int ptn){
	int TableId = TbInlet_p[PipeId];
	float T = Table_step[TableId][ptn];
	for(int i=0;i<pts_p[PipeId];i++){
		int pointId = pId_p[PipeId][i];
		tem_pp[pointId][ptn] = T;
	}
}

void pipe_film(float* Q_Right,int e,int pn,int ptn,int pointId){
	int num = NodeNum_e[e];
	float T1 = 0, T2 = 0;

	T2 = tem_pp[pointId][ptn];

	fprintf(log_bound,"\n水管边界条件 pn = %d\n",pn);
	for(int k=0;k<num;k++){
		if(pn==0){
			if (k<8)	T1=t_e[e];
			else		T1=0.0;
		}else{
			if (k<8)	T1 = tem_pp[pointId][ptn-1];
			else		T1 = tem_pp[pointId][ptn-1];
		}
		Q_Right[k] -= (T1+T2)*Q3[e][k]*beta_p;
		fprintf(log_bound,"Q_Right [%d] : %f\n",k,Q_Right[k]);
	}
}