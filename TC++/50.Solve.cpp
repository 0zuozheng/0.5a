/* * * * * * * * * * * * * * * * * *
*
*		求解线性方程组
*
* * * * * * * * * * * * * * * * * */

#include "SOLVE.h"

//定义计算变量
//int nlca,nuca;
int nonzero;
float *B, *X;
extern int conn;
extern void NodWrite(int,int);
extern int *NumCols, **ColId;
extern float **ValueRowCol;

void Tstep(int ln,int ptn,int pn, float dt) //Called by Displacement.  ln -> loadcase_now;  ptn -> period_totalnow; pn -> period_now; dt -> dalt
{ 
	int i,j,k,row,col,m;
	float begintime,endtime,value,beginvalue,endvalue,tmp, *Qtotal;	

	B=new float[equationnum];
	for(i=0;i<equationnum;i++)  B[i] = 0.0;
	for(i=0;i<equationnum;i++)	NumCols[i] = 0;
	nonzero=0;
	for(i=0;i<elementnum_c;i++){
		int num = class_e[i];
		Qtotal = (float *)calloc(num,sizeof(float));
		m=material_e[i];
		begintime=time_e[i];  endtime=begintime+dt;
		for(j=0;j<num;j++) Qtotal[j]=0;
		if(deactive_l[i]==false){
			for(j=0;j<num;j++){
				if(node_e[i][j]>=0){
					row=node_e[i][j]-order[node_e[i][j]];
					for(k=0;k<num;k++){
						if(node_e[i][k]>=0){
							col=node_e[i][k]-order[node_e[i][k]];
							value=H[i][j][k]+2*P[i][j][k]/dt;
							Insert(row,col,value);
						}
					}
				}
			}
			
			for(j=0;j<6;j++){ // 第三类边界条件
				if (surf_e[i][j]>0){
					if(SurfLoadType[(surf_e[i][j])]==3)	third_film(Qtotal,i,j,ln,pn,ptn,surf_e[i][j]);
				}
			}
			if(begintime<200){  beginvalue=exp(-m_m[m]*begintime);  endvalue=exp(-m_m[m]*endtime); }
			else beginvalue=endvalue=0;

			if(m_m[m]!=0)	// 水化热
				for(j=0;j<num;j++)
					Qtotal[j]-=(beginvalue+endvalue)*Q0[i][j];	//组合Q的水化热项（上一步+本步）,注意这里的符号为负
			
			//else ???
			for(j=0;j<num;j++){
				tmp=0;
				for(k=0;k<num;k++){
					if(deactive_l[i]==false && time_e[i]==0)
						 tmp += (H[i][j][k]-2*P[i][j][k]/dt)*t_e[i];
					else tmp += (H[i][j][k]-2*P[i][j][k]/dt)*t_n[node_e[i][k]];
				}	Qtotal[j] += tmp;
			}

			for (j=0;j<num;j++)	fprintf(DEBUG,"%d - %d : %d - %f\n",ln,ptn+1,j+1,Qtotal[j]);

			for(j=0;j<num;j++){
				if(calculate_n[node_e[i][j]]==true){
					col=node_e[i][j]-order[node_e[i][j]];
					B[col]-=Qtotal[j];							//B放在等号右边，符号为负
				}
			}
		}free(Qtotal);
	}

	//for(i=0;i<boundarynum_c;i++){	// 消行修正法 《传热学的有限元方法》 P67
	//	for(j=0;j<num_b[i];j++){
	//		int nod = list_b[i][j];
	//		if(calculate_n[nod]==true){
	//			col = nod-order[nod];
	//			printf("!!! fix id = %d\n",col);
	//			for(int k=0;k<nodenum_c;k++){
	//				ValueRowCol[col][k] = 0;
	//				ValueRowCol[k][col] = 0;
	//			}	ValueRowCol[col][col] = 1;

	//			B[col] = fix_b[i];
	//		}
	//	}
	//}

	for(i=0;i<boundarynum_c;i++){	// 对角线扩大法 《传热学的有限元方法》 P67
		for(j=0;j<num_b[i];j++){
			int nod = list_b[i][j];
			if(calculate_n[nod]==true){
				col = nod-order[nod];
				for (int k=0;k<=col;k++){
					if (ColId[col][k]==col){
						ValueRowCol[col][k] *= 1e10f;
					}
				}
				B[col]=ValueRowCol[col][col]*fix_b[i];
			}
		}
	}

	//for(i=0;i<fixnum_c;i++){
	//	if(node_b[i]<maxnode_l[ln] && calculate_n[node_b[i]]==true){
	//		col=node_b[i]-order[node_b[i]];
	//		AD[col]=AD[col]*1e10f;  
	//        if (abs(fix_b[i])>1e-6) B[col]=AD[col]*fix_b[i];
	//	}
	//}

	fprintf(log_solver,"step = %d, period = %d, dalt = %f\n",ln,pn+1,dt);
	SolveWrite(ln,pn);
	Tecplot(ln,pn,ptn);		//Tecplot
	NodWrite(ptn,2);
}

void Insert(int row, int col, float value)
{
	if(row>=col){	// 位于下三角阵中
		int NewCol=0;
		for(int i=0;i<NumCols[row];i++){
			if(ColId[row][i]==col){
				ValueRowCol[row][i] += value;
				NewCol = 1;
			}
		}if(NewCol==0){
			int count = NumCols[row];
			nonzero++;
			if(count>conn)	printf("半带宽大于要求！\n");
			NumCols[row] += 1;
			ColId[row][count] = col;
			ValueRowCol[row][count] = value;
		}
	}
}

void SolveWrite(int ln,int pn) //ln -> loadcase_now; pn -> period_now
{
	int i,j,k=0;

	Imsl_f_sparse_elem *A;

	A=new Imsl_f_sparse_elem[nonzero];

	fprintf(log_solver,"\n--------A---------\n");
	for(i=0;i<equationnum;i++){
		fprintf(log_solver,"%8d:",i);
		for(j=0;j<NumCols[i];j++){
			A[k].row=i;
			A[k].col=ColId[i][j];
			A[k].val=ValueRowCol[i][j];
			if (A[k].val>1000)	fprintf(log_solver,"%8d %8.1e,",A[k].col,A[k].val);
			else	fprintf(log_solver,"%8d%8.5f,",A[k].col,A[k].val);
			k++;
		}fprintf(log_solver,"\n");
	}

	fprintf(log_solver,"--------B--------\n");
	for (i=0;i<equationnum;i++)
		fprintf(log_solver,"%2d%16.5f\n",i+1,B[i]);
	fprintf(log_solver,"\n----------X--------\n");

	X = imsl_f_lin_sol_posdef_coordinate (equationnum, k, A, B, 0);

	for (i=0;i<equationnum;i++)
		fprintf(log_solver,"%2d%16.5f\n",i+1,X[i]);

	for(i=0;i<nodenum_c;i++)
		if (marked[i]>0) t_n[i]=X[i-order[i]];
		else t_n[i]=999.99f;
	delete[] B;		delete[] A;
	fprintf(log_answer,"Temperature: loadcase=%d,period=%d\n",ln,pn+1);
	for(i=0;i<nodenum_c;i++){
		if(t_n[i]==999.99f)	fprintf(log_answer,"%4d %16s\n",i+1,"NULL");
		else				fprintf(log_answer,"%4d %16.3f\n",i+1,t_n[i]);
	}
	if(rst_flag)	Rst_wrt_node(nodenum_c,t_n);	//写出NODES结果
}