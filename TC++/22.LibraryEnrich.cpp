#include "ELEMENT.h"

extern float Rb_p;
extern float Modulus(float*);
extern void Normalize(float*);
extern int PointPerTri[8];
extern FILE *log_enrich;

float LMN_Node8[10][3] = // LMN:代表局部坐标
{{-1,-1,-1},
 { 1,-1,-1},
 { 1, 1,-1},
 {-1, 1,-1},
 {-1,-1, 1},
 { 1,-1, 1},
 { 1, 1, 1},
 {-1, 1, 1}};

int NodeOnFace[6][4] = 
{{3,2,1,0},
 {4,5,6,7},
 {4,0,1,5},
 {5,1,2,6},
 {6,2,3,7},
 {7,3,0,4}};

int Prism35[4][6] =
{{0,8,1,3,9,2},
 {4,8,0,7,9,3},
 {1,8,5,2,9,6},
 {5,8,4,6,9,7}};

int Prism46[4][6] = 
{{1,8,2,0,9,3},
 {5,8,1,4,9,0},
 {2,8,6,3,9,7},
 {6,8,5,7,9,4}};

int Prism12[4][6] =
{{3,8,2,7,9,6},
 {2,8,1,6,9,5},
 {1,8,0,5,9,4},
 {0,8,3,4,9,7}};

int Axis_Of_Face[66];
int PointPrism[6];
int PointPrismX4[4][6];
float LMN_Prism[6][3];

void write_element_Hammer(int en, int in, int classnum,float r,float dir[3]){	// in: integral point
	fprintf(log_N,"\n*********  Gauss Point [%d]  *********\n",in);
	float xx=0,yy=0,zz=0;
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

void GetLMNofPrism(int Prism[4][6],int no,int *NodeList){
	for (int j=0;j<6;j++){
		PointPrism[j] = NodeList[Prism[no][j]];
		PointPrismX4[no][j] = PointPrism[j];
		for(int k=0;k<3;k++){
			LMN_Prism[j][k]=LMN_Node8[Prism[no][j]][k];
		}
	}
}

void GetPrismLMN(int CrossFace,int no,int *NodeList){	// Triangular prism: 2 triangles and 3 squares
	int No = 0;
	if (CrossFace==35)	GetLMNofPrism(Prism35,no,NodeList);
	if (CrossFace==46)	GetLMNofPrism(Prism46,no,NodeList);
	if (CrossFace==12)	GetLMNofPrism(Prism12,no,NodeList);
}

void circular(float XYZ[3],float P1[3],float P2[3],float *r,float d[3]){ // 一定要在同一坐标系下，即水管结点和积分点都在全局坐标系
	float PI = 3.1415926;	float XYZ_0[3]; float sum=0.0;
	float gx = XYZ[0]-P1[0],  gy = XYZ[1]-P1[1], gz = XYZ[2]-P1[2];
	float lx = P2[0]-P1[0],   ly = P2[1]-P1[1],  lz = P2[2]-P1[2];
	float A = (lx*gx+ly*gy+lz*gz)/(lx*lx+ly*ly+lz*lz);
	XYZ_0[0] = A*lx + P1[0];
	XYZ_0[1] = A*ly + P1[1];
	XYZ_0[2] = A*lz + P1[2];
	for(int i=0;i<3;i++)
		d[i] = XYZ[i] - XYZ_0[i];
	*r = Modulus(d);
	Normalize(d);
	*r /= Rb_p;
}

void Enrichsf(float x,float y,float z,float r,float *sf_){ // 丰富形函数
	GetSf(x,y,z,sf_);	// 前8项

	float b = 6.4f;	// ZUOZUOZ
	float phi = - exp(-b*r); // f(r)=-exp(-b*r), b=6.4
	
	for (int i=0;i<8;i++)	sf_[i+8] = phi*sf_[i];	// 扩充形函数
}

void Enrichdxyzdsf(float x,float y, float z, float r,float *dir,float *sf_,float **dsf_,float *det_,float **dxydsf_){
	float **dxydsf_t, dphi[3];
	Alloc2DArray_float(&dxydsf_t,3,8);
	GetDetDxyzsf(det_,dsf_,dxydsf_t);

	float b=6.4;
	float phi = - exp(-b*r);

	extern float Rb_p;

	dphi[0] = b*exp(-b*r)*dir[0]/Rb_p;
	dphi[1] = b*exp(-b*r)*dir[1]/Rb_p;
	dphi[2] = b*exp(-b*r)*dir[2]/Rb_p;

	

	for (int i=0;i<3;i++){
		for(int j=0;j<8;j++){
			dxydsf_[i][j] = dxydsf_t[i][j];
			dxydsf_[i][j+8] = phi*dxydsf_[i][j] + dphi[i]*sf_[j];
		}
	}
}

int N_Hammer(int en, int order, int GaussNum, float Qi[][4], float Wi[]){	// 管点坐标已经含在node_e[8],node_e[9]中
	int CrossFace = 12;
	int TriNum = 4;
	int PointNum = TriNum * GaussNum;

	det[en] = (float *)calloc(PointNum,sizeof(float));		w[en] = (float *)calloc(PointNum,sizeof(float));
	Alloc2DArray_float(&sf[en],PointNum,8*2);
	Alloc3DArray_float(&dsf[en],PointNum,3,8);				Alloc3DArray_float(&dxyzsf[en],PointNum,3,8*2);

	LMN_Pipe_e[en][0][0] = 0;	LMN_Pipe_e[en][0][1] =  0;	LMN_Pipe_e[en][0][2] = -1;
	LMN_Pipe_e[en][1][0] = 0;	LMN_Pipe_e[en][1][1] =  0;	LMN_Pipe_e[en][1][2] =  1;
	
	for(int i=0;i<2;i++){	// 两个管结点
		for(int j=0;j<3;j++){
			LMN_Node8[i+8][j] = LMN_Pipe_e[en][i][j];
		}
	}
	
	GetEleXYZ(10,en);	// 获取结点坐标
	for(int i=0;i<10;i++){
		fprintf(log_enrich,"%d, %f, %f, %f\n", node_e[en][i]+1, nodeXYZ_e[i][0], nodeXYZ_e[i][1], nodeXYZ_e[i][2]);
	}

	for(int i=0;i<TriNum;i++){
		GetPrismLMN(CrossFace,i,node_e[en]);
		for(int j=0;j<6;j++){
			for(int k=0;k<3;k++){
				fprintf(DEBUG,"%f\t",LMN_Prism[j][k]);
			}	fprintf(DEBUG,"\n");
		}

		float volume = 2.0 / 2;

		for(int j=0;j<GaussNum;j++){
			int PointID = i*GaussNum + j;
			float r,dir[3],LMN_local[2][3]={0.0},XYZ_global[3]={0.0},coor[3];

			for(int k=0; k<2; k++)			//先在两个三角形上确定积分点位置
				for(int p=k*3; p<(k+1)*3; p++) //三个点
					for(int x=0;x<3;x++) //三个坐标分量
						LMN_local[k][x] += Qi[j][p%3] * LMN_Prism[p][x];

			for(int l=0;l<3;l++)
				coor[l] = (1 - Qi[j][3])*LMN_local[0][l]/2 + (1 + Qi[j][3])*LMN_local[1][l]/2;	// 线性插值

			GetXYZ(coor,nodeXYZ_e,XYZ_global);		// 获得全局坐标
			circular(XYZ_global,nodeXYZ_e[8],nodeXYZ_e[9],&r,dir); // 计算全局坐标下的半径
						
			Enrichsf(coor[0],coor[1],coor[2],r,sf[en][PointID]);
			GetDsf(coor[0],coor[1],coor[2],dsf[en][PointID]);
			Enrichdxyzdsf(coor[0],coor[1],coor[2],r,dir,sf[en][PointID],dsf[en][PointID],&det[en][PointID],dxyzsf[en][PointID]);
			w[en][PointID] = volume * Wi[j];

			write_element_Hammer(en,PointID,16,r,dir);
			fprintf(log_enrich,"**, %f, %f, %f\n",coor[0],coor[1],coor[2]);	// 输出积分点 -> 结点
			fprintf(log_enrich,"%d, %f, %f, %f\n",100+PointID,XYZ_global[0],XYZ_global[1],XYZ_global[2]);	// 输出积分点 -> 结点
		}		
	}
	fprintf(log_enrich,"*Element, type=C3D6\n");
	for(int i=0;i<4;i++){
		fprintf(log_enrich,"%d,",i+1);
		for(int j=0;j<6;j++){
			fprintf(log_enrich,"%d,",PointPrismX4[i][j]+1);
		}	fprintf(log_enrich,"\n");
	}fclose(log_enrich);
	log_enrich = fopen("enrich.inp","a");
	return PointNum;	// 斜穿情况：分成包含顶、末点的一个大六面体 + 三个小六面体组成的L型体，共剖分5段求解
}

void ShapeFunc_enrich(int en){
	PointNum_e[en] = plan_e[en] * PointPerTri[plan_e[en]];
	if (plan_e[en]==2)	PointNum_e[en] = N_Hammer(en,2,PointNum_e[en],ti2xi2,tw2xw2);
	if (plan_e[en]==3)	PointNum_e[en] = N_Hammer(en,3,PointNum_e[en],ti3xi3,tw3xw3);
	if (plan_e[en]==4)	PointNum_e[en] = N_Hammer(en,4,PointNum_e[en],ti4xi4,tw4xw4);
	if (plan_e[en]==5)	PointNum_e[en] = N_Hammer(en,5,PointNum_e[en],ti5xi5,tw5xw5);
}