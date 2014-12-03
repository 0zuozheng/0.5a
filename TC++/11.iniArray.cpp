/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *   Array/Matrix Initialization Tool
 *
 *	 Author:          Z Zuo
 *   Last Modified :  Dec. 2014
 *
 * (All the funcs can be called in Data.h)
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdlib.h>

/* Calloc a float matrix named arr with k1*k2 zeroes ... */
void Alloc2DArray_float(float*** arr,int k1,int k2){
	*arr = (float **)calloc(k1,sizeof(float *));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (float *)calloc(k2,sizeof(float));
	}
}

/* Realloc a float matrix named arr from k0*k2 to k1*k2 ... */
void Realloc2DArray_float(float ***arr, int k0, int k1, int k2){
	*arr = (float **)realloc(*arr, (k1)*sizeof(float*));
	for (int i=k0;i<k1;i++){
		(*arr)[i] = (float *)calloc(k2,sizeof(float));
	}
}

/* Calloc an integer matrix named arr with k1*k2 zeroes ... */
void Alloc2DArray_int (int*** arr,int k1,int k2){
	*arr = (int **)calloc(k1,sizeof(int *));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (int *)calloc(k2,sizeof(int));
	}
}

/* Calloc a string matrix named arr with k1*k2 spaces ... */
void Alloc2DArray_char (char*** arr,int k1,int k2){
	*arr = (char **)calloc(k1,sizeof(char *));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (char *)calloc(k2,sizeof(char));
	}
}

/* Calloc a float matrix named arr with k1*k2*k3 zeroes ... */
void Alloc3DArray_float(float**** arr,int k1,int k2,int k3){
	*arr = (float ***)calloc(k1,sizeof(float **));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (float **)calloc(k2,sizeof(float *));
		for (int j=0;j<k2;j++){
			(*arr)[i][j] = (float *)calloc(k3,sizeof(float));
		}
	}
}

/* Calloc an integer matrix named arr with k1*k2*k3 zeroes ... */
void Alloc3DArray_int  (int**** arr,int k1,int k2,int k3){
	*arr = (int ***)calloc(k1,sizeof(int **));
	for (int i=0;i<k1;i++){
		(*arr)[i] = (int **)calloc(k2,sizeof(int *));
		for (int j=0;j<k2;j++){
			(*arr)[i][j] = (int *)calloc(k3,sizeof(int));
		}
	}
}

/* Free an ineger matrix named arr with num rows ... */
void FREE2D_int(int** arr,int num){
	for(int i=0;i<num;i++)
		free(arr[i]);
	free(arr);
}

/* Free a float matrix named arr with num rows ... */
void FREE2D_float(float** arr,int num){
	for(int i=0;i<num;i++)
		free(arr[i]);
	free(arr);
}

/* Free a string matrix named arr with num strings ... */
void FREE2D_char(char** arr,int num){
	for(int i=0;i<num;i++)
		free(arr[i]);
	free(arr);
}

/* Free an ineger matrix named arr with num1*num2 rows ... */
void FREE3D_int(int*** arr,int num1,int num2){
	for(int i=0;i<num1;i++){
		for(int j=0;j<num2;j++){
			free(arr[i][j]);
		}free(arr[i]);
	}free(arr);
}

/* Free a float matrix named arr with num1*num2 rows ... */
void FREE3D_float(float*** arr,int num1,int num2){
	for(int i=0;i<num1;i++){
		for(int j=0;j<num2;j++){
			free(arr[i][j]);
		}free(arr[i]);
	}free(arr);
}