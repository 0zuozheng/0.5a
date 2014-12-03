/* * * * * * * * * * * * * * * * * *
*
*	 刚度矩阵安装
*
* * * * * * * * * * * * * * * * * */

#include "CONSTRUCT.h"

//定义计算变量
float *time_e, **Q0, ***H, ***P, ***PP, **Q3;
float time_totalnow = 0.0;
int *NumCols, **ColId;
float **ValueRowCol;
int conn=120;		// !!! 120 !!!
int *elements_n, equationnum, *marked, *order;
bool *calculate_n; 
extern FILE *log_answer,*log_check,*out,*inp_step;
extern int SurfLoadId;
ofstream tecplot;	// 结果输出文件流

extern void Read_step(FILE*);
extern void restart_files();
extern void close_files();
extern void FreeVariables();
extern void FreeShape();

void FreeAssembly(){
	free(time_e);		free(calculate_n);		free(elements_n);
	FREE2D_float(Q0,elementnum_c);			FREE3D_float(H,elementnum_c,8);
	FREE3D_float(P,elementnum_c,8);			FREE2D_float(Q3,elementnum_c);
	
	FREE2D_int(ColId,enrich_nodenum_c);		FREE2D_float(ValueRowCol,enrich_nodenum_c);
	free(NumCols);
}

void Info_step(int step){
	for(int j=0;j<40;j++)	fprintf(log_check,"-");
	fprintf(log_check,"\n#### Step %d \n",step);
	fprintf(log_check,"\n Step Name : %s\n Total time : %f\n Total increments : %d\n\n",name_l,totaltime_l,inc_l);
	fprintf(log_check,"Deactivated element list :\n");
	for(int i=0,j=0;j<elementnum_c;j++){
		if ((i%8==0)&&(i!=0))	fprintf(log_check,"\n");
		if (deactive_l[j]){
			fprintf(log_check,"%8d",j+1);
			i++;
		}
	}
	fprintf(log_check,"\nFace information :\n");
	for(int i=0;i<elementnum_c;i++){
		if ((i%8==0)&&(i!=0))	fprintf(log_check,"\n");
		for(int j=0;j<6;j++){
			fprintf(log_check,"%d",surf_e[i][j]);
		}fprintf(log_check,"  ");
	}fprintf(log_check,"\n");

	fprintf(log_check,"\nInteractive Conditions :\n");
	for(int i=1;i<SurfLoadId;i++){
		fprintf(log_check,"  Id = %d Type = %d beta = %f TableID = %d\n",i,SurfLoadType[i],BetaValue[i],SinkTable[i]+1);
	}fprintf(log_check,"\n");
	
	fprintf(log_check,"\nBoundary Conditions :\n");
	for(int i=0;i<boundarynum_c;i++){
		fprintf(log_check,"  Id = %d Number of Nodes = %d\n",i+1,num_b[i]);
		if(table_b[i]!=-1)	fprintf(log_check,"  Boundary condition is depent on the table No.%d\n",table_b[i]+1);
		for(int j=0;j<num_b[i];j++){
			fprintf(log_check,"\t%d",list_b[i][j]+1);
		}
	}
}

float TimeInter(float time,int i,int j,int num,float* timetable,float* valuetable){
	for (int i=0;i<num-1;i++){
		if ((timetable[i]<=time)&&(timetable[i+1]>=time)){
			float d = (time-timetable[i])/(timetable[i+1]-timetable[i]);
			float f = valuetable[i] + d*(valuetable[i+1]-valuetable[i]);
			return f;
		}
	}
}

void AdjustTable(float time_now,int ptn){
	for (int i=0;i<tablenum_c;i++){
		Table_step[i] = (float*)realloc(Table_step[i],(ptn+inc_l)*sizeof(float));
		for (int j=0;j<inc_l;j++){
			float time = time_now + timeinc_l * (j+0.5f);	//线性假定
			Table_step[i][ptn+j] = TimeInter(time,i,j,num_tb[i],time_tb[i],value_tb[i]);
		}
	}fprintf(DEBUG,"#### Table for increment=%d time=%f \n",ptn,time_now);
	for(int i=0;i<tablenum_c;i++)
		for(int j=0;j<inc_l+ptn;j++)
			fprintf(DEBUG,"Inc = %d Value = %f \n",j,Table_step[i][j]);

	for (int i=0;i<ppnum_c;i++)	// tem_pp插值
		tem_pp[i] = (float*)realloc(tem_pp[i],(ptn+inc_l)*sizeof(float));
}

void AssignP(int pn){
	for(int i=0;i<elementnum_c;i++){
		if(pipe_e[i]){
			for (int j=0;j<16;j++){
				for (int k=0;k<16;k++){
					P[i][j][k] = PP[i][j][k];
				}
			}
			if (pn==0){
				for(int k=0;k<8;k++)	for(int l=8;l<16;l++)	P[i][k][l]=0;
				for(int k=8;k<16;k++)	for(int l=0;l<8;l++)	P[i][k][l]=0;
			}
		}
		else{
			for (int j=0;j<8;j++)
				for (int k=0;k<8;k++)
					P[i][j][k] = PP[i][j][k];
		}
	}
}

void Assembly(){		// Called By Main.cpp
	char ansname[100];	tecplot.open("tec.dat");
	int step = 0, period_totalnow = 0;
	maxelement_l = elementnum_c;

	calculate_n = (bool*)calloc(enrich_nodenum_c,sizeof(bool));	elements_n = (int*)calloc(enrich_nodenum_c,sizeof(int));
	t_n = (float*)calloc(enrich_nodenum_c,sizeof(float));
	time_e = (float*)calloc(elementnum_c,sizeof(int));			deactive_l	= (int*)calloc(elementnum_c,sizeof(int));
	P  = (float***)calloc(elementnum_c,sizeof(float**));		PP = (float***)calloc(elementnum_c,sizeof(float**));
	Q0 = (float**)calloc(elementnum_c,sizeof(float*));			H = (float***)calloc(elementnum_c,sizeof(float**));
	Q3 = (float**)calloc(elementnum_c,sizeof(float*));

	NumCols     = (int*)calloc(enrich_nodenum_c,sizeof(int));
	Alloc2DArray_int(&ColId,enrich_nodenum_c,conn);				Alloc2DArray_float(&ValueRowCol,enrich_nodenum_c,conn);
	Alloc2DArray_int(&surf_e,elementnum_c,6);					Alloc2DArray_float(&Table_step,tablenum_c,1);

	printf("Now Starting Calculation...\n");

	while(!feof(inp_step)){
		step ++;	restart_files();
		Read_step(inp_step);
		sprintf(ansname,"RESULT %03d(%s).dat",step,name_l);
		log_answer = fopen(ansname,"w");
		marked = (int*)calloc(enrich_nodenum_c,sizeof(int));		order  = (int*)calloc(enrich_nodenum_c,sizeof(int));
		Info_step(step);	// 输出信息
		AdjustTable(time_totalnow,period_totalnow);		// 表格插值
		PreTreatment(step);
		for(int j=0;j<inc_l;j++){
			cout<<"  Step = "<<step<<" Increment = "<<j+1<<" under calculating......"<<endl;
			time_totalnow += timeinc_l;	// 注意：步长一致
			
			AssignP(j);
			Tstep(step,period_totalnow,j,timeinc_l);

			period_totalnow ++;
			for(int k=0;k<elementnum_c;k++){
				if(deactive_l[k] == false)	time_e[k] += timeinc_l;
			}
		}fclose(log_answer);
		free(marked);    free(order);
	}close_files();
	FreeShape();	// 释放形函数
	FreeVariables();// 释放全局变量
	FreeAssembly();	// 释放线性方程组矩阵
}