#include "Solve.h"
void NodWrite(int ptn, int point){
	char cmd[500];
	sprintf(cmd,"Node=%d.rpt",point+1);
	FILE *fp = fopen(cmd,"a");
	if (ptn==0)	fprintf(fp,"X = %f, Y = %f, Z = %f\n",xyz_n[point][0],xyz_n[point][1],xyz_n[point][2]);
	fprintf(fp,"%8d%10.4f\n",ptn,t_n[point]);
	fclose(fp);
}