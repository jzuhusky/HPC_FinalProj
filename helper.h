/*
 * Helper Function File
 */

#include "util.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

int powerOfTwo(int N);
void printMat(int N, double **u);
void printMatPlot(int N, double **u, int rank, int procsPerRow);
void zeroOut(int N, double **u);
void addMatrices(double **A, double **B, int dim);
double** myMalloc(int N);
void myFree(double **mat, int N);

int powerOfTwo(int N){
	if ( N && (N & (N-1)) == 0 ){
		return 1;
	}else{ 
		return 0;
	}
}

void printMat(int N, double **u){
	int i,j;
	printf("***********PRINT MAT************\n\n");
	for(i=0;i<N;i++){
		printf("[");
		for(j=0;j<N;j++){
			printf(" %.9f ", u[i][j]);
		}
		printf(" ]\n");
	}
}

void printMatPlot(int N, double **u, int rank, int procsPerRow){
	int i,j;
	double h = 1.0/(N*procsPerRow-1.0*procsPerRow);
	if (rank == 0 ){
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				//if ( u[i][j] == 0.0 ){printf("Rank %d (i,j) = (%d,%d)\n",rank,i,j);}
				//printf("%f %f %f\n",h*((rank/procsPerRow)*N+i),h*((rank%procsPerRow)*N+j),u[i][j]);
				printf("%02f %02f %f #%d\n",h*((rank/procsPerRow)*N+i),h*((rank%procsPerRow)*N+j),u[i][j],rank);
			}
		}
	}else if (rank % procsPerRow == 0 ){
		for(i=1;i<N;i++){
			for(j=0;j<N;j++){
				//if ( u[i][j] == 0.0 ){printf("Rank %d (i,j) = (%d,%d)\n",rank,i,j);}
				//printf("%f %f %f\n",h*((rank/procsPerRow)*N+i),h*((rank%procsPerRow)*N+j),u[i][j]);
				printf("%02f %02f %f #%d\n",h*((rank/procsPerRow)*N+(i-rank/procsPerRow)),h*((rank%procsPerRow)*N+j),u[i][j],rank);
			}
		}
	}else if ( rank < procsPerRow  && rank > 0 ){
		for(i=0;i<N;i++){
			for(j=1;j<N;j++){
				//if ( u[i][j] == 0.0 ){printf("Rank %d (i,j) = (%d,%d)\n",rank,i,j);}
				//printf("%f %f %f\n",h*((rank/procsPerRow)*N+i),h*((rank%procsPerRow)*N+j),u[i][j]);
				printf("%02f %02f %f #%d\n",h*((rank/procsPerRow)*N+i),h*((rank%procsPerRow)*N+(j-rank%procsPerRow)),u[i][j],rank);
			}
		}
	}else{
		for(i=1;i<N;i++){
			for(j=1;j<N;j++){
				//if ( u[i][j] == 0.0 ){printf("Rank %d (i,j) = (%d,%d)\n",rank,i,j);}
				//printf("%f %f %f\n",h*((rank/procsPerRow)*N+i),h*((rank%procsPerRow)*N+j),u[i][j]);
				printf("%02f %02f %f #%d\n",h*((rank/procsPerRow)*N+(i-rank/procsPerRow)),h*((rank%procsPerRow)*N+(j-rank%procsPerRow)),u[i][j],rank);

			}
		}
	}

}

void zeroOut(int N, double **u){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			u[i][j] = 0.0;
		}
	}
}


double** myMalloc(int N){
	int i;
	double **mat = (double**)calloc(N,sizeof(double*));
	for(i=0;i<N;i++){
		mat[i] = (double*)calloc(N,sizeof(double));
	}
	return mat;
}

// Helper free function for 2D matricies....
void myFree(double **mat, int N){
	int i;
	for(i=0;i<N;i++){
		free(mat[i]);
	}
	free(mat);
}

void addMatrices(double **A, double **B, int dim){
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			A[i][j] = A[i][j]+B[i][j];
		}
	}
}
