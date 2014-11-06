/* * * * * * * * * * * * * * * * * *
*
*	 刚度矩阵安装
*
* * * * * * * * * * * * * * * * * */

#include "CONSTRUCT.h"

//定义计算变量

float *time_e, **Q0, ***H, ***P, **Q3;
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
	
	FREE2D_int(ColId,nodenum_c);			FREE2D_float(ValueRowCol,nodenum_c);
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
			float time = time_now + timeinc_l * (j+0.5);	//线性假定
			Table_step[i][ptn+j] = TimeInter(time,i,j,num_tb[i],time_tb[i],value_tb[i]);
		}
	}fprintf(DEBUG,"#### Table for increment=%d time=%f \n",ptn,time_now);
	for(int i=0;i<tablenum_c;i++){
		for(int j=0;j<inc_l+ptn;j++){
			fprintf(DEBUG,"Inc = %d Value = %f \n",j,Table_step[i][j]);
		}
	}
}

void Assembly(){		// Called By Main.cpp
	tecplot.open("tec.dat");
	char ansname[100];
	int step = 0, period_totalnow = 0;
	maxelement_l = elementnum_c;

	time_e = (float*)calloc(elementnum_c,sizeof(int));		deactive_l	= (int*)calloc(elementnum_c,sizeof(int));
	P  = (float***)calloc(elementnum_c,sizeof(float**));	//PP = (float***)calloc(elementnum_c,sizeof(float**));
	Q0 = (float**)calloc(elementnum_c,sizeof(float*));		H = (float***)calloc(elementnum_c,sizeof(float**));
	Q3 = (float**)calloc(elementnum_c,sizeof(float*));	
	calculate_n = (bool*)calloc(nodenum_c,sizeof(bool));	elements_n=(int*)calloc(nodenum_c,sizeof(int));

	NumCols     = (int*)calloc(nodenum_c,sizeof(int));
	Alloc2DArray_int(&ColId,nodenum_c,conn);
	Alloc2DArray_float(&ValueRowCol,nodenum_c,conn);

	Alloc2DArray_int(&surf_e,elementnum_c,6);
	Alloc2DArray_float(&Table_step,tablenum_c,1);

	printf("Now Starting Calculation...\n");

	while(!feof(inp_step)){
		step ++;	restart_files();
		Read_step(inp_step);
		sprintf(ansname,"RESULT %03d(%s).dat",step,name_l);
		log_answer = fopen(ansname,"w");
		marked = (int*)calloc(nodenum_c,sizeof(int));
		order  = (int*)calloc(nodenum_c,sizeof(int));
		Info_step(step);	// 输出信息
		AdjustTable(time_totalnow,period_totalnow);		// 表格插值
		PreTreatment(step);
		for(int j=0;j<inc_l;j++){
			cout<<"  Step = "<<step<<" Increment = "<<j+1<<" under calculating......"<<endl;
			time_totalnow += timeinc_l;	// 注意：步长一致
			Tstep(step,period_totalnow,j,timeinc_l);

			period_totalnow++;
			for(int k=0;k<elementnum_c;k++){
				if(deactive_l[k]==false)	time_e[k]+=timeinc_l;
			}
		}fclose(log_answer);
		free(marked);    free(order);
	}close_files();
	FreeShape();	// 释放形函数
	FreeVariables();// 释放全局变量
	FreeAssembly();	// 释放线性方程组矩阵
}

void PrintElement(int ln,int i){
	int ii,jj;
	fprintf(log_element,"\nLoadcase = %d, Element = %d\n",ln,i+1);
	fprintf(log_element,"--------P -----------\n");
	for (ii=0;ii<8;ii++){
		fprintf(log_element,"%2d:",ii+1);
		for (jj=0;jj<8;jj++)	fprintf(log_element,"%8.3f",P[i][ii][jj]);
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
}

void PreTreatment(int ln){ //ln -> loadcase_now

	for(int i=0;i<nodenum_c;i++) { marked[i]=1;  elements_n[i]=0;  calculate_n[i]=false;}  //marked 标记该结点激活
	for(int i=0;i<elementnum_c;i++){		
		if(!deactive_l[i]){	// active -> 单元激活
			for(int j=0;j<8;j++){
				elements_n[node_e[i][j]]++;	calculate_n[node_e[i][j]]=true;
			}
			for(int j=0;j<6;j++){
				if (surf_e[i][j]>0)	Tface(ln,i,j,surf_e[i][j]);	// 面组装
			}
			
			StiffMatrix_HEX8(i);		// 形成单元刚度矩阵

			PrintElement(ln,i);
		}
	}
	equationnum=nodenum_c;		int l=0;
	for(int i=0;i<nodenum_c;i++){
		if(calculate_n[i]==false) {	equationnum--;	l++;  marked[i]=0; }
		order[i]=l;
	}
}