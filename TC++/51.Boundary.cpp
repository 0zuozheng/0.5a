#include "Construct.h"

void third_film(float *Q2,int e,int f,int ln,int pn,int ptn,int LoadId){ // e -> element; f -> face; ln -> loadcase_now; pn -> period_now; ptn -> period_totalnow;
	float T1=0, T2=0;
	int   tableid = SinkTable[LoadId];
	float beta    = BetaValue[LoadId];

	T2 = Table_step[tableid][ptn];

	if(deactive_l[e]==false && pn==0)
		T1 = t_e[e];
	else
		T1 = Table_step[tableid][ptn-1];

	for(int k=0;k<8;k++){
		Q2[k] -= (T1+T2)*Q3[e][k]*beta;
	}

	fprintf(DEBUG,"%d - %f\n",ptn,T2);
}