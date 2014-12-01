#include "ELEMENT.h"
extern FILE *log_enrich;
extern int PointPrismX4[4][6];
extern float LMN_Prism[6][3];
extern int PointPerTri[8];
extern float LMN_Node8[10][3];
extern void Enrichsf(float[3],float,float*);
extern void Enrichdxyzdsf(float,float*,float*,float**,float*,float**);
extern void circular(float[3],float[3],float[3],float*,float[3]);
extern int GetCrossDir(float**);
extern void GetPrismLMN(int,int,int*);

void write_element_Hammer(int en, int in, int classnum, float r, float dir[3]){	// in: integral point
	fprintf(log_N,"\n*********  Hammer Point [%d]  *********\n",in);
	float xx=0, yy=0, zz=0;
	for (int k=0;k<8;k++){
		xx += sf[en][in][k]*nodeXYZ_e[k][0];	yy += sf[en][in][k]*nodeXYZ_e[k][1];	zz += sf[en][in][k]*nodeXYZ_e[k][2];
	}fprintf(log_N,"Location：x = %7.4f, y = %7.4f, z = %7.4f\n",xx,yy,zz);
	fprintf(log_N,"\tr = %7.4f dir = [%7.4f %7.4f %7.4f]\n",r,dir[0],dir[1],dir[2]);

	fprintf(log_N,"ShapeFunc sf [1~%d]:\n\t",classnum);
	for(int k=0;k<classnum;k++)	fprintf(log_N,"%7.4f\t",sf[en][in][k]);
	fprintf(log_N,"\nDerivative: dsf");
	for (int k=0;k<3;k++){
		fprintf(log_N,"\n\t");
		for (int l=0;l<classnum/2;l++)	fprintf(log_N,"%7.4f\t",dsf[en][in][k][l]);		
	}fprintf(log_N,"\nDerivative: dxydsf");
	for (int k=0;k<3;k++){
		fprintf(log_N,"\n\t");
		for (int l=0;l<classnum;l++)	fprintf(log_N,"%7.4f\t",dxyzsf[en][in][k][l]);		
	}
	fprintf(log_N,"\nWeight w : %7.4f\n",w[en][in]);
	fprintf(log_N,"Det det : %7.4f\n",det[en][in]);	
}

void AddPipeNode(int en,float **LMN){		float XYZ[3];
	extern FILE *out;	fprintf(out,"*Node\n");
	int n = pipe_nodenum_c;
	Realloc2DArray_float(&xyz_n, n, n+2, 3);	
	for (int i=0;i<2;i++){
		node_e[en][i+8] = n+i;			fprintf(out,"%d ",n+i+1);
		GetXYZ(LMN[i],nodeXYZ_e,XYZ);
		for(int j=0;j<3;j++){
			xyz_n[n+i][j] = XYZ[j];		fprintf(out,", %.2f ",XYZ[j]);
		}fprintf(out,"\n");
	}
	pipe_nodenum_c += 2;

	for(int i=0;i<2;i++)	// 2个管结点
		for(int j=0;j<3;j++)// 3个方向
			LMN_Node8[i+8][j] = LMN[i][j];
}

void Write_PipeTopo(int en){
/*** 记录水管拓扑到enrich.inp中等待Abq读入验证 **/	
	fprintf(log_enrich,"*Node\n");
	for(int i=0;i<10;i++)	fprintf(log_enrich,"%d, %f, %f, %f\n", node_e[en][i]+1, nodeXYZ_e[i][0], nodeXYZ_e[i][1], nodeXYZ_e[i][2]);	
	fprintf(log_enrich,"*Element, type=C3D6\n");	
	for(int i=0;i<4;i++){
		fprintf(log_enrich,"%d,",i+1);
		for(int j=0;j<6;j++){
			fprintf(log_enrich,"%d,",PointPrismX4[i][j]+1);
		}	fprintf(log_enrich,"\n");
	}fclose(log_enrich);
	log_enrich = fopen("enrich.inp","a");
}

void Get_r_dir(int PointID, float Q[4], float coor[3], float *r, float dir[3]){
	float LMN_local[2][3]={0.0},XYZ_global[3]={0.0};
	for(int k=0; k<2; k++)			//先在两个三角形上确定积分点位置
		for(int p=k*3; p<(k+1)*3; p++) //三个点
			for(int x=0;x<3;x++) //三个坐标分量
				LMN_local[k][x] += Q[p%3] * LMN_Prism[p][x];
	for(int l=0;l<3;l++)
		coor[l] = (1 - Q[3])*LMN_local[0][l]/2 + (1 + Q[3])*LMN_local[1][l]/2;	// 线性插值

	GetXYZ(coor,nodeXYZ_e,XYZ_global);		// 获得全局坐标
	circular(XYZ_global, nodeXYZ_e[8], nodeXYZ_e[9], r, dir); // 计算全局坐标下的半径
	fprintf(log_enrich,"%d, %f, %f, %f\n",10000+PointID,XYZ_global[0],XYZ_global[1],XYZ_global[2]);	// 输出积分点 -> 结点
}

int N_Hammer(int en, int order, int GaussNum, float Qi[][4], float Wi[]){	// 管点坐标已经含在 LMN_Pipe_e[en][0/1][0/1/2] 中, Qi包含三个面积坐标，一个轴向坐标
	int TriNum = 4;		int PointNum = TriNum * GaussNum;
	det[en] = (float *)calloc(PointNum,sizeof(float));		w[en] = (float *)calloc(PointNum,sizeof(float));
	Alloc2DArray_float(&sf[en],PointNum,8*2);
	Alloc3DArray_float(&dsf[en],PointNum,3,8);				Alloc3DArray_float(&dxyzsf[en],PointNum,3,8*2);

	GetEleXYZ(8,en);	// 获取结点坐标
	LMN_Pipe_e[en][0][0] = 0;	LMN_Pipe_e[en][0][1] =  0;	LMN_Pipe_e[en][0][2] = -1;	// ZUOZUO 这里的先后顺序以后要有逻辑化
	LMN_Pipe_e[en][1][0] = 0;	LMN_Pipe_e[en][1][1] =  0;	LMN_Pipe_e[en][1][2] =  1;	// ZUOZUO 通过这个数组实际上就可以判断Cross方向
	int CrossFace = GetCrossDir(LMN_Pipe_e[en]);

	AddPipeNode(en,LMN_Pipe_e[en]);	// ZUOZUO 这步在未来取消，本来这个信息就应该在PipeGeometry中解决
	GetEleXYZ(10,en);

	for(int i=0;i<TriNum;i++){
		float volume = 2.0 / 2;	// ZUOZUO 这步要换成体积计算
		GetPrismLMN(CrossFace,i,node_e[en]);
		for(int j=0;j<GaussNum;j++){
			float r, dir[3], LMN[3];
			int PointID = i*GaussNum + j;
			Get_r_dir(PointID, Qi[j], LMN, &r, dir);
			Enrichsf(LMN,r,sf[en][PointID]);
			GetDsf(LMN[0],LMN[1],LMN[2],dsf[en][PointID]);
			Enrichdxyzdsf(r,dir,sf[en][PointID],dsf[en][PointID],&det[en][PointID],dxyzsf[en][PointID]);
			w[en][PointID] = volume * Wi[j];
			write_element_Hammer(en,PointID,16,r,dir);
		}
	}Write_PipeTopo(en);

	for (int j=0;j<8;j++){ // 最后把这个单元的结点都标记， 注意：前面的node_e里面[8][9]是水管结点，现在换成富集结点
		if (!enrichorder_n[node_e[en][j]])	// 标记富集结点
				enrichorder_n[node_e[en][j]] = (enrich_nodenum_c++);	// 结点编号从0开始
		node_e[en][j+8] = enrichorder_n[node_e[en][j]];
	}
	return PointNum;	// 斜穿情况：分成包含顶、末点的一个大六面体 + 三个小六面体组成的L型体，共剖分5段求解
}

void ShapeFunc_enrich(int en){
	PointNum_e[en] = plan_e[en] * PointPerTri[plan_e[en]];
	if (plan_e[en]==2)	PointNum_e[en] = N_Hammer(en,2,PointNum_e[en],ti2xi2,tw2xw2);
	if (plan_e[en]==3)	PointNum_e[en] = N_Hammer(en,3,PointNum_e[en],ti3xi3,tw3xw3);
	if (plan_e[en]==4)	PointNum_e[en] = N_Hammer(en,4,PointNum_e[en],ti4xi4,tw4xw4);
	if (plan_e[en]==5)	PointNum_e[en] = N_Hammer(en,5,PointNum_e[en],ti5xi5,tw5xw5);
}