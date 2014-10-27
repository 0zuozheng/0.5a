#include <stdlib.h>

void  Alloc2DArray_float(float*** arr,int k1,int k2){
	*arr = (float **)calloc(k1,sizeof(float *));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (float *)calloc(k2,sizeof(float));
	}
}

void    Alloc2DArray_int  (int*** arr,int k1,int k2){
	*arr = (int **)calloc(k1,sizeof(int *));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (int *)calloc(k2,sizeof(int));
	}
}

void    Alloc2DArray_char  (char*** arr,int k1,int k2){
	*arr = (char **)calloc(k1,sizeof(char *));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (char *)calloc(k2,sizeof(char));
	}
}

void	Alloc3DArray_float(float**** arr,int k1,int k2,int k3){
	*arr = (float ***)calloc(k1,sizeof(float **));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (float **)calloc(k2,sizeof(float *));
		for (int j=0;j<k2;j++){
			(*arr)[i][j] = (float *)calloc(k3,sizeof(float));
		}
	}
}

void   Alloc3DArray_int  (int**** arr,int k1,int k2,int k3){
	*arr = (int ***)calloc(k1,sizeof(int **));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (int **)calloc(k2,sizeof(int *));
		for (int j=0;j<k2;j++){
			(*arr)[i][j] = (int *)calloc(k3,sizeof(int));
		}
	}
}

void FREE2D_int(int** arr,int num){
	for(int i=0;i<num;i++)
		free(arr[i]);
	free(arr);
}

void FREE2D_float(float** arr,int num){
	for(int i=0;i<num;i++)
		free(arr[i]);
	free(arr);
}

void FREE2D_char(char** arr,int num){
	for(int i=0;i<num;i++)
		free(arr[i]);
	free(arr);
}

void FREE3D_int(int*** arr,int num1,int num2){
	for(int i=0;i<num1;i++){
		for(int j=0;j<num2;j++){
			free(arr[i][j]);
		}free(arr[i]);
	}free(arr);
}

void FREE3D_float(float*** arr,int num1,int num2){
	for(int i=0;i<num1;i++){
		for(int j=0;j<num2;j++){
			free(arr[i][j]);
		}free(arr[i]);
	}free(arr);
}