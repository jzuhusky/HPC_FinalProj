/*
 *	Joe Zuhusky
 *	Serial Implementation of 2D Multigrid-Method
 *	for 2D Poission Equation
 *	First Attempt at A simple serial multigrid 
 *	Implementation
 */

#include "util.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define v1 3
#define v2 3

int relaxCount=0;

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

void printMatPlot(int N, double **u){
	int i,j;
	double h = 1.0/(N+1.0);
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf("%d %d %f\n",i,j,u[i][j]);
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

void interpolateToFine(int fromDim, double **coarseGrid, double **fineGrid);
void restrictToCoarse(int fromDim,double **dest, double **origin);
void computeResidual(int dim, double **u, double **res, double **rho);
double residualNorm(int dim, double **u, double **res, double **rho);
void relax(int iterations, double **u,double **rhs, int dim);
void addMatrices(double **A, double **B, int dim);
double** myMalloc(int N);
void myFree(double **mat, int N);

int main(int argc, char **argv){
	if ( argc < 3 ){
		printf("Please enter the correct number of parameters\n");
		exit(1);
	}

	int numVcycles, i, j, k, N, numberOfGrids, n2,C;
	N             = atoi(argv[1]); // Input Initial Grid Size (Use power of two)

	if ( argc == 3 ){
		numberOfGrids = atoi(argv[2]);
	}else{
		numberOfGrids = floor(log2(N));
	}
	if ( powerOfTwo(N) == 0 ){
		printf("Please use a Dim that is a power of 2\n");
		exit(0);
	}else if ( pow(2,numberOfGrids) > N ){
		printf("Please use a smaller number of grids, #grids <= log2(N)\n");
		exit(0);
	}
	N   = N + 1;
	// Data Arrays for Pointers to Grids
	double **u, **defect, **defect2;
	double **temp, **uPtrs[numberOfGrids], **defects[numberOfGrids], **rhs[numberOfGrids], **tempDefect;

	// Preprocessing -> Allocate Space for all of the Grids/Defects
	n2 = N;
	for(i=0;i<numberOfGrids;i++){
		uPtrs[i]   = myMalloc(n2);
		rhs[i]     = myMalloc(n2);
		defects[i] = myMalloc(n2);
		zeroOut(n2,uPtrs[i]);
		zeroOut(n2,defects[i]);
		zeroOut(n2,rhs[i]);
		n2 = n2 / 2 + 1;
	}
	// Initialize Source Function
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			rhs[0][i][j] = 1.0; // Simple Source Function...	
		}
	}

	// Preprocessing -> Restrict Right hand Side(Source Function) to coarser Grids
	// Again, just allocating some more memory for these grids
	double **current;
	double r_norm, initNorm;
	double conv_tol = 0.00000001;
	// For each C, do a V-Cycle
	timestamp_type time1,time2;
	get_timestamp(&time1);
	r_norm = residualNorm(N,uPtrs[0],defects[0], rhs[0]);
	initNorm=r_norm;
	printf("%f\n",r_norm / initNorm);
	while ( (r_norm / initNorm) > conv_tol ){
		r_norm = residualNorm(N,uPtrs[0],defects[0], rhs[0]);
		printf("Norm=%0.8f -> %8f %% of init\n",r_norm, r_norm/initNorm * 100);
		for (k=0; k<numberOfGrids-1;k++){
			relax(v1,uPtrs[k],rhs[k], N);
			computeResidual(N,uPtrs[k],defects[k],rhs[k]);
			restrictToCoarse(N,rhs[k+1],defects[k]);
			N = N / 2 + 1;
			zeroOut(N,uPtrs[k+1]);
		}
		relax(50*v2,uPtrs[numberOfGrids-1],rhs[numberOfGrids-1], N);
		for (k=numberOfGrids-1; k > 0;k--){
			current = myMalloc(2*N-1); // Some Temp Memory to holder interpolated value...
                        interpolateToFine(N,uPtrs[k],current);
			N = 2*N-1;
			addMatrices(uPtrs[k-1],current,N);
			relax(v2,uPtrs[k-1],rhs[k-1],N);
			myFree(current,N); // Find a way to get rid of this...
		}
	} // END Full Multigrid Loop
	//printf("Relaxation count: %d\n",relaxCount);
//	printMatPlot(N,uPtrs[0]);
	/*Clean Up and Finish */
	// TODO FREE ALL MEMORY USED!!!!
	get_timestamp(&time2);
        double elapsed = timestamp_diff_in_seconds(time1,time2);
	printf("# Grids: %d\n",numberOfGrids);
        printf("Total Time %f seconds.\n", elapsed);

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

void interpolateToFine(int fromDim, double **coarseGrid, double **fineGrid){
	int i,j;
	int N = fromDim*2 - 1;
	int state;
	// Using Bi-linear Interpolation
	for(i=1;i<N-1;i++){
		for(j=1;j<N-1;j++){
			if ( i%2 == 0 && j%2 == 0 ) state = 0; 		
			else if ( i%2 == 1 && j%2 == 0 ) state = 1; 		
			else if ( i%2 == 0 && j%2 == 1 ) state = 2; 		
			else if ( i%2 == 1 && j%2 == 1 ) state = 3;
			switch(state){
			case (0): // These Points are direct transforms from Coarse to Fine
				fineGrid[i][j] = coarseGrid[i/2][j/2];
				break;
			case (1): 
				// Odd Row, Even Columns
				// Average Points Up and Down
				fineGrid[i][j] = (coarseGrid[i/2][j/2]+coarseGrid[i/2+1][j/2])/2.0;
				break;
			case (2):
				// Even Row, Odd Columns
				// Same Idea as Case 1
				fineGrid[i][j] = (coarseGrid[i/2][j/2]+coarseGrid[i/2][j/2+1])/2.0;
				break;
			case (3):
				// Odd Row Odd Columns
				// Average 4 Coarse Points around this point...
				fineGrid[i][j] = (coarseGrid[i/2][j/2]+coarseGrid[i/2+1][j/2]+coarseGrid[i/2][j/2+1]+coarseGrid[i/2+1][j/2+1])/4.0;
				break;
			} 		
		}
	}
}

void restrictToCoarse(int fromDim,double **dest, double **origin){
	int i,j;
	int N = fromDim/2 + 1;
	// options for Resctrion 
	// 1. Injection
	// 2. Full Weighing
	// 3. Half Weighing
	// USING HERE: Half!! Weighing
	// First take care of boundaries...

	for(i=0;i<N;i++){
		dest[i][0] = origin[2*i][0];
		dest[i][N-1] = origin[2*i][N-1];
		dest[0][i] = origin[0][2*i];
		dest[N-1][i] = origin[N-1][i];
	}

	for(i=1;i<N-1;i++){
		for(j=1;j<N-1;j++){
			dest[i][j] = 1.0/8.0* (
				origin[2*i][2*j+1]+origin[2*i][2*j-1]+origin[2*i+1][2*j]+origin[2*i-1][2*j]+4.0*origin[2*i][2*j]
			); 
		}
	}
}


void computeResidual(int dim, double **u, double **res, double **rho){
	int i,j, N=dim;
	double h = 1.0/(N + 1.0);
	for(i=1;i<dim-1;i++){
	        for(j=1; j<dim-1; j++){
			res[i][j] = (rho[i][j] - (4.0*u[i][j]-u[i+1][j]-u[i-1][j]-u[i][j+1]-u[i][j-1])/h/h);
		}
	}
}

double residualNorm(int dim, double **u, double **res, double **rho){
	int i,j, N=dim;
	double h = 1.0/(N + 1.0);

	double resD=0, temp=0;

	for(i=1;i<dim-1;i++){
	        for(j=1; j<dim-1; j++){
			temp = (rho[i][j] - (4.0*u[i][j]-u[i+1][j]-u[i-1][j]-u[i][j+1]-u[i][j-1])/h/h); 
			resD += temp * temp;
		//	resD += res[i][j]*res[i][j];
		}
	}
	return sqrt(resD);

}

void relax(int iterations, double **u,double **rhs, int dim){
	relaxCount++; // Macro Counter
	int i,j, iter, procsPerRow, N=dim;
	double h = 1.0/(N+1.0);
	double w = 1.0;
	for(iter=0;iter<iterations;iter++){
		for (i=1; i<dim-1; i++){
			for (j=1; j < dim-1; j++){
				u[i][j] = u[i][j] + w * (u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1]+rhs[i][j]*h*h-4.0*u[i][j])/4.0;
			}
		}
		for (i=dim-2; i> 0; i--){
			for (j=dim-2; j > 0; j--){
				u[i][j] = u[i][j] + w * (u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1]+rhs[i][j]*h*h-4.0*u[i][j])/4.0;
			}
		}
	}
}
