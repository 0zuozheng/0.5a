#include "CONSTRUCT.h"

extern void LineSource(int, int);

void PrintElement(int ln,int i){
	int ii,jj;
	fprintf(log_element,"\nLoadcase = %d, Element = %d\n",ln,i+1);
	fprintf(log_element,"--------PP -----------\n");
	for (ii=0;ii<8;ii++){
		fprintf(log_element,"%2d:",ii+1);
		for (jj=0;jj<8;jj++)	fprintf(log_element,"%8.3f",PP[i][ii][jj]);
		fprintf(log_element,"\n");
	};
	fprintf(log_element,"--------Q0 -----------\n");
	for (ii=0;ii<8;ii++)
		fprintf(log_element,"%2d:%8.3f\n",ii+1,Q0[i][ii]);
	fprintf(log_element,"--------H -----------\n");
	for (ii=0;ii<8;ii++) {
		fprintf(log_element,"%2d",ii+1);
		for (jj=0;jj<8;jj++) fprintf(log_element,"%8.3f",H[i][ii][jj]);
		fprintf(log_element,"\n");
	};
	fprintf(log_element,"--------Q3 -----------\n");
	for (ii=0;ii<8;ii++) 
		fprintf(log_element,"%2d:%8.3f\n",ii+1,Q3[i][ii]);
	
	fclose(log_element);
	log_element = fopen("Element.log","a");
}

void PreTreatment(int ln){ //ln -> loadcase_now
	for(int i=0;i<enrich_nodenum_c;i++) { marked[i]=1;  elements_n[i]=0;  calculate_n[i]=false;}  //marked 标记该结点激活
	for(int i=0;i<elementnum_c;i++){
		int num;
		if(!deactive_l[i]){	// active -> 单元激活
			if(pipe_e[i])	num = 16;
			else			num = 8;
			for(int j=0;j<num;j++){
				elements_n[node_e[i][j]]++;	calculate_n[node_e[i][j]]=true;
			}
			for(int j=0;j<6;j++){
				if (surf_e[i][j]>0)	Tface(ln,i,j,surf_e[i][j]);	// 面组装
			}
			if(pipe_e[i])	StiffMatrix_HEX16(i);		// 形成单元刚度矩阵
			else			StiffMatrix_HEX8(i);		// 形成单元刚度矩阵

			if(pipe_e[i])	LineSource(i, 2);	// 管壁边界条件

			PrintElement(ln,i);
		}
	}
	equationnum = enrich_nodenum_c;		int l=0;
	for(int i=0;i<enrich_nodenum_c;i++){
		if(calculate_n[i]==false) {	equationnum--;	l++;  marked[i]=0; }
		order[i]=l;
	}
}