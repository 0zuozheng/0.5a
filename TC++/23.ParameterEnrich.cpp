#include "ELEMENT.h"

int PointPrism[6];
int PointPrismX4[4][6];
float LMN_Prism[6][3];

float LMN_Node8[10][3] = // LMN:代表局部坐标
{{-1,-1,-1},
 { 1,-1,-1},
 { 1, 1,-1},
 {-1, 1,-1},
 {-1,-1, 1},
 { 1,-1, 1},
 { 1, 1, 1},
 {-1, 1, 1},
 { 0, 0, 0},
 { 0, 0, 0}};

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

void GetLMNofPrism(int Prism[4][6],int no,int *NodeList){
	for (int j=0;j<6;j++){
		PointPrism[j] = NodeList[Prism[no][j]];
		PointPrismX4[no][j] = PointPrism[j];
		for(int k=0;k<3;k++){
			LMN_Prism[j][k]=LMN_Node8[Prism[no][j]][k];
		}
	}
}

void GetPrismLMN(int CrossAxis,int no,int *NodeList){	// Triangular prism: 2 triangles and 3 squares
	int No = 0;	
	if (CrossAxis==0)	GetLMNofPrism(Prism46,no,NodeList);	// X axis Face 4,6
	if (CrossAxis==1)	GetLMNofPrism(Prism35,no,NodeList);	// Y axis Face 3,5
	if (CrossAxis==2)	GetLMNofPrism(Prism12,no,NodeList); // Z axis Face 1,2
}

int GetCrossDir(float **LMN){
	if (LMN[0][0]*LMN[1][0]==-1)	return 0;
	else if (LMN[0][1]*LMN[1][1]==-1)	return 1;
	else if (LMN[0][2]*LMN[1][2]==-1)	return 2;
	else return -1;
}