#include "ELEMENT.h"

extern float Modulus(float*);
extern void Normalize(float*);

float beta_p = 53663.7f;		// lmd = 1, beta = lmd/(1.2*r0*ln(r0/rp)), 185/(1.2*0.84/100*log(1.2))
float Rb_p   = 0.84f;	// 对应1.5m×1.5m
float x_p = 0.65f, y_p = 0.45f;	// 中心

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
float LMN_Prism[6][3];

void write_element_Hammer(int en, int in, int classnum){	// in: integral point
	fprintf(log_N,"\n\tGauss Point [%d]\n",in);
	float xx=0,yy=0,zz=0;
	for (int k=0;k<8;k++){
		xx += sf[en][in][k]*nodeXYZ_e[k][0];	yy += sf[en][in][k]*nodeXYZ_e[k][1];	zz += sf[en][in][k]*nodeXYZ_e[k][2];
	}fprintf(log_N,"\t\tLocation：x = %f, y = %f, z = %f\n",xx,yy,zz);

	fprintf(log_N,"\t\tShapeFunc sf [1~%d]:\n",classnum);
	for(int k=0;k<classnum;k++)	fprintf(log_N,"%f\t",sf[en][in][k]);
	fprintf(log_N,"\n\t\tDerivative: dsf\n");
	for (int k=0;k<3;k++){
		for (int l=0;l<classnum/2;l++)	fprintf(log_N,"%f\t",dsf[en][in][k][l]);
		fprintf(log_N,"\n");
	}fprintf(log_N,"\t\tDerivative: dxydsf\n");
	for (int k=0;k<3;k++){
		for (int l=0;l<classnum;l++)	fprintf(log_N,"%f\t",dxyzsf[en][in][k][l]);
		fprintf(log_N,"\n");
	}
}

void GetLMNofPrism(int Prism[4][6],int no,int *NodeList){
	for (int j=0;j<6;j++){
		PointPrism[j] = NodeList[Prism[no][j]];
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

void GetXYZ(float coor[3],float (*XYZ_e)[3],float XYZ[3]){
	float sf_[8];
	GetSf(coor[0],coor[1],coor[2],sf_);
	for(int i=0;i<3;i++)	XYZ[i] = 0.0f;

	for (int i=0;i<8;i++)	XYZ[0] += sf_[i]*XYZ_e[0][i];
	for (int i=0;i<8;i++)	XYZ[1] += sf_[i]*XYZ_e[1][i];
	for (int i=0;i<8;i++)	XYZ[2] += sf_[i]*XYZ_e[2][i];
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
	fprintf(DEBUG,"det = %f\n",*det_);

	float b=6.4;
	float phi = - exp(-b*r);

	extern float Rb_p;

	dphi[0] = b*exp(-b*r)*dir[0]/Rb_p;
	dphi[1] = b*exp(-b*r)*dir[1]/Rb_p;
	dphi[2] = b*exp(-b*r)*dir[2]/Rb_p;

	fprintf(log_N,"r = %f dir = [%f %f %f] dphi = [%f %f %f]\n",r,dir[0],dir[1],dir[2],dphi[0],dphi[1],dphi[2]);

	for (int i=0;i<3;i++){
		for(int j=0;j<8;j++){
			dxydsf_[i][j] = dxydsf_t[i][j];
			dxydsf_[i][j+8] = phi*dxydsf_[i][j] + dphi[i]*sf_[j];
		}
	}
}

void N_Hammer(int en, int order, int GaussNum, float Qi[][3], int CrossFace){	// 管点坐标已经含在node_e[8],node_e[9]中
	order = 3; CrossFace = 12;	GaussNum = 3*4;
	int axis = Axis_Of_Face[CrossFace];
	int TriNum = 4;	int PointNum = TriNum * GaussNum;

	det[en] = (float *)calloc(PointNum,sizeof(float));		Alloc2DArray_float(&sf[en],PointNum,8*2);
	Alloc3DArray_float(&dsf[en],PointNum,3,8);				Alloc3DArray_float(&dxyzsf[en],PointNum,3,8*2);

	float i_tri[100][3],i_line[100];
	LMN_Node8[8][0] = 0;	LMN_Node8[8][1] =  0;	LMN_Node8[8][2] = -1;
	LMN_Node8[9][0] = 0;	LMN_Node8[9][1] =  0;	LMN_Node8[9][2] =  1;

	GetEleXYZ(10,en);	// 获取结点坐标
	for(int i=0;i<3;i++)
		for(int j=0;j<8;j++)
			fprintf(DEBUG,"%d %d = %f\n",i,j,nodeXYZ_e[j][i]);

	PrismOrder(order, i_tri, i_line);

	for(int i=0;i<TriNum;i++){
		GetPrismLMN(CrossFace,i,node_e[en]);

		for(int j=0;j<GaussNum;j++){
			int PointID = i*GaussNum+j;
			float r,dir[3],LMN_local[2][3]={0.0},XYZ_global[3]={0.0},coor[3];			

			for(int k=0; k<2; k++)			//先在两个三角形上确定积分点位置
				for(int p=k*3; p<(k+1)*3; p++) //三个点
					for(int x=0;x<3;x++) //三个坐标分量
						LMN_local[k][x] += i_tri[j][p%3] * LMN_Prism[p][x];

			for(int l=0;l<3;l++)
				coor[l] = (1 - i_line[j])*LMN_local[0][l]/2 + (1 + i_line[j])*LMN_local[1][l]/2;	// 线性插值

			GetXYZ(coor,nodeXYZ_e,XYZ_global);		// 获得全局坐标
			circular(XYZ_global,nodeXYZ_e[8],nodeXYZ_e[9],&r,dir);// 计算全局坐标下的半径
						
			Enrichsf(coor[0],coor[1],coor[2],r,sf[en][PointID]);
			GetDsf(coor[0],coor[1],coor[2],dsf[en][PointID]);
			Enrichdxyzdsf(coor[0],coor[1],coor[2],r,dir,sf[en][PointID],dsf[en][PointID],&det[en][PointID],dxyzsf[en][PointID]);

			write_element_Hammer(en,j,16);

			fprintf(DEBUG,"Node8 = %f, %f, %f\n",nodeXYZ_e[8][0],nodeXYZ_e[8][1],nodeXYZ_e[8][2]);
			fprintf(DEBUG,"Node9 = %f, %f, %f\n",nodeXYZ_e[9][0],nodeXYZ_e[9][1],nodeXYZ_e[9][2]);
			fprintf(DEBUG,"r = %f, dir = %f, %f, %f\n",r,dir[0],dir[1],dir[2]);
			fprintf(DEBUG,"%d, %f, %f, %f\n",100+PointID,XYZ_global[0],XYZ_global[1],XYZ_global[2]);
		}
		fprintf(DEBUG,"%d,",i+1);
		for (int j=0;j<6;j++)	fprintf(DEBUG,"%d,",PointPrism[j]+1);
		fprintf(DEBUG,"\n");
	}
	// 斜穿情况：分成包含顶、末点的一个大六面体 + 三个小六面体组成的L型体，共剖分5段求解
}

void ShapeFunc_enrich(int en){
	if (plan_e[en]==1008)	N_Hammer(en,0, 8,i2x2x2_,35);
	if (plan_e[en]==1027)	N_Hammer(en,0,27,i3x3x3_,35);
	if (plan_e[en]==1064)	N_Hammer(en,0,64,i4x4x4_,35);
	if (plan_e[en]==1125)	N_Hammer(en,0,125,i5x5x5_,35);
}