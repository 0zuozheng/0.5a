//#include	"Data.h"
//#include	"Transform.h"
//extern FILE* log_transform;
//
//int node_face[6][4] = { {1,2,6,5},
//						{0,4,7,3},
//						{2,3,7,6},
//						{0,1,5,4},
//						{0,3,2,1},
//						{4,5,6,7}};
//
//void SinglePipeCrossFace(){
//
//	int k1,k2,k3;
//	float Xaxis[3] = {1, 0, 0};	// 新坐标系中向量（X轴）
//	float B[4][3] = {0};		// 存储ABCD在新坐标系下的坐标
//	float coor_new[3];			// 交汇点原坐标系内坐标	
//	float pipe[3],XofO;			// 存储管线方向
//	float N[3];					// 存储插值函数
//	
//	float norm = NomalizeVector(pipe,L1,L2);		// 计算管线向量
//	fprintf(log_transform,"Norm = %f\n",norm);
//
//	CalcRotate(pipe,Xaxis);	// 旋转坐标系，使水管与新坐标系下的X轴平行
//	
//	CalcTranslate(L1);		// 平移坐标系，使水管位于原点，即水管变为X轴
//	
//	for (int i=0;i<4;i++){	// 计算A平移后的四点坐标，并录入到B
//		CalcCoordinate(A[i],B[i],rotatinMatrix,translateMatrix);
//	}
//	int N1 = JudgePointInQuad(B[0],B[1],B[2],B[3]);	// 判断点是否在四边形内
//	
//	if (N1 >  0){
//		float  N3 = JudgePointInTri(A[0],A[1],A[2]);
//		if (N3>0)	k1 = 0, k2 = 1, k3 = 2;	
//		else		k1 = 2, k2 = 3, k3 = 0;	
//		
//		TriInterFunc(N,A[k1],A[k2],A[k3]);
//		float XofABC[3] = {B[k1][0],B[k2][0],B[k3][0]};
//		XofO = Interpolate(XofABC, N, 3);
//		if ((XofO<=norm)&&(XofO>=0)){
//			fprintf(log_transform,"交汇点 x=%f 在管线内部\n",XofO);
//			CoorInterpolateByTri(coor_new,A[k1],A[k2],A[k3],N);
//		}else	fprintf(log_transform,"交汇点 x=%f 在管线外部\n",XofO);
//	}
//
//	if (N1 == 0){
//		// ToDoToDoToDoToDoToDoToDoToDoToDo
//		//    在边上的有待进一步编程处理
//		// ToDoToDoToDoToDoToDoToDoToDoToDo
//		fprintf(log_transform,"点在四边形边上");
//	}
//	if (N1 <  0) fprintf(log_transform,"点在四边形外");
//}
//
//void PipeCrossFace(int PipeId,int EleId,int FaceId){
//	
//	
//	
//	for (int i=0;i<pipenum_c;i++){
//		for (int j=0;j<pts_p[i]-1;j++){
//			SinglePipeCrossFace(xyz_p[i][j],xyz_p[i][j+1]);
//		}
//	}
//}
//
//void GeometryAnalysis(){
//	for (int k=0;k<pipenum_c;k++){
//		for (int i=0;i<elementnum_c;i++){
//			for (int j=0;j<6;j++){
//				PipeCrossFace(k,i,j);
//			}
//		}
//	}
//}
//
