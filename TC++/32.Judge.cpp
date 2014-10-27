#include "Transform.h"

int JudgePointInQuad(float* A,float* B,float* C,float* D){	// O为原点，不再给出;ABCD逆时针排序
	float N1 = 1;	float* ABCD[4];
	ABCD[0] = A;	ABCD[1] = B;	ABCD[2] = C;	ABCD[3] = D;

	for (int i=0;i<4;i++){
		N1 *= PlaneCrossProduct(ABCD[i],ABCD[(i+1)%4],0);	// k = 0, 代表投影向X轴
		printf("Product[%d] = %f\n",i,N1);
		if (N1==0){
			if((ABCD[i][1]==ABCD[(i+1)%4][1])&&(ABCD[i][2]==ABCD[(i+1)%4][2])){
				N1=1;				
			}else{
				N1 = JudgePointInLine(ABCD[i],ABCD[(i+1)%4]);	
				break;
			}
		}
	}
	printf("N1=%f",N1);

	if (N1>0){		fprintf(log_transform,"\n点在四边形内\n");	return  1;	}
	else if (N1==0){fprintf(log_transform,"\n点在四边形上\n");	return  0;	}
	else{			fprintf(log_transform,"\n点在四边形外\n");	return -1;	}
}

float JudgePointInLine(float* A,float* B){	// 默认A[0]=B[0]=0.0
	float N2;

	N2 = A[1]*B[1] + A[2]*B[2];	//	N2<0, 代表在线段上; N2=0, 代表点在端点上; N2>0, 代表点在线段外
	return (-N2-abs(N2));		//  N2<=0, 返回0; N2>0, 返回负数
}

float JudgePointInTri(float* A,float* B,float* C){	// O为原点，不再给出;ABCD逆时针排序
	float N3 = 1;

	N3 = (C[2]-A[2])/(C[1]-A[1])-(B[2]-A[2])/(B[1]-A[1]);
	N3 *= (C[2]-A[2])/(C[1]-A[1])-A[2]/A[1];

	if (N3>0)	fprintf(log_transform,"\n点在三角形ABC内\n");
	else		fprintf(log_transform,"\n点在三角形ABC外\n");

	return N3;
}