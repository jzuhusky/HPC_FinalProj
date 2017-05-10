/*
 *	Joe Zuhusky
 *	Serial Implementation of 3D Multigrid-Method
 *	for 3D Poission Equation
 *	Implementation
 *
 *
 *	WARNING: THIS CODE IS NOT YET WORKING!!!!!!!!!!!
 *
 *
 */

#include "util.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define v1 3
#define v2 3

int relaxCount=0;

void printMatPlot(int N, double ***mat){

	int i,j,k;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf("%d %d %f\n",i,j,mat[i][j][N/2]);
		}
	}

}

int powerOfTwo(int N){
	if ( N && (N & (N-1)) == 0 ){
		return 1;
	}else{ 
		return 0;
	}
}

void zeroOut(int N, double ***u){
	int i,j,k;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(j=0;j<N;j++){
				u[i][j][k] = 0.0;
			}
		}
	}
}
double*** myMalloc(int N);
void addMatrices(double ***A, double ***B, int dim);
void myFree(double ***mat, int N);
void relax(int iterations, double ***u,double ***rhs, int dim);
double residualNorm(int dim, double ***u, double ***res, double ***rho);
void computeResidual(int dim, double ***u, double ***res, double ***rho);
void interpolateToFine(int fromDim, double ***coarseGrid, double ***fineGrid);
void restrictToCoarse(int fromDim,double ***dest, double ***origin);

int main(int argc, char **argv){
	if ( argc != 3 ){
		printf("Please enter the correct number of parameters\nParam#1->Dimension Size (Power of 2), Param#2 = Num V-Cyclesn\n");
		exit(1);
	}
	int i, j, k, N, numberOfGrids, n2,C;
	N             = atoi(argv[1]); // Input Initial Grid Size (Use power of two)
	numberOfGrids = floor(log2(N));
	if ( powerOfTwo(N) == 0 ){
		printf("Please use a Dim that is a power of 2\n");
		exit(0);
	}else if ( pow(2,numberOfGrids) > N ){
		printf("Please use a smaller number of grids, #grids <= log2(N)\n");
		exit(0);
	}
	N   = N + 1;
	// Data Arrays for Pointers to Grids
	double ***uPtrs[numberOfGrids], ***defects[numberOfGrids], ***rhs[numberOfGrids];

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
	// Initialize BCs
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				rhs[0][i][j][k] = 0.0;
			}
		}
	}
	rhs[0][N/2][N/2][N/2] = 1.0;
	rhs[0][N/2][N/2][N/2+1] = 1.0;
	rhs[0][N/2][N/2+1][N/2] = 1.0;
	rhs[0][N/2][N/2+1][N/2+1] = 1.0;
	rhs[0][N/2+1][N/2][N/2] = 1.0;
	rhs[0][N/2+1][N/2][N/2+1] = 1.0;
	rhs[0][N/2+1][N/2+1][N/2] = 1.0;
	rhs[0][N/2+1][N/2+1][N/2+1] = 1.0;

	// Preprocessing -> Restrict Right hand Side(Source Function) to coarser Grids
	// Again, just allocating some more memory for these grids
	double ***current;
	double r_norm, initNorm;
	double conv_tol = 0.000001;
	// For each C, do a V-Cycle
	printf("N = %d\n",N);
	r_norm = residualNorm(N,uPtrs[0],defects[0], rhs[0]);
	initNorm=r_norm;
	for(C=0; C<atoi(argv[2]); C++){
		r_norm = residualNorm(N,uPtrs[0],defects[0], rhs[0]);
		printf("Norm: %8f %% of initial\n",r_norm/initNorm * 100);
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
			relax(v2,uPtrs[k-1],rhs[k-1], N);
			myFree(current,N); // Find a way to get rid of this...
		}
	} // END Full Multigrid Loop
	//printf("Relaxation count: %d\n",relaxCount);
	//printMatPlot(N,uPtrs[0]);
}

double*** myMalloc(int N){
	int i,j;
	double ***mat = (double***)calloc(N,sizeof(double**));
	for(i=0;i<N;i++){
		mat[i] = (double**)calloc(N,sizeof(double*));
		for (j=0;j<N;j++){
			mat[i][j] = (double*)calloc(N,sizeof(double));
		}
	}
	return mat;
}

// Helper free function for 3D matricies....
void myFree(double ***mat, int N){
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			free(mat[i][j]);
		}
		free(mat[i]);
	}
	free(mat);
}

void addMatrices(double ***A, double ***B, int dim){
	int i,j,k;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			for(k=0;k<dim;k++){
				A[i][j][k] += B[i][j][k];
			}
		}
	}
}

void interpolateToFine(int fromDim, double ***coarseGrid, double ***fineGrid){
	int i,j,k;
	int N = fromDim*2 - 1;
	int state;
	// Using Tri-linear Interpolation
	for(i=1;i<N-1;i++){
		for(j=1;j<N-1;j++){
			for(k=1;k<N-1;k++){	
				if ( i%2 == 0 && j%2 == 0 && k%2==0) state = 0; 		
				else if ( i%2 == 1 && j%2 == 0 && k%2 == 0 ) state = 1; 		
				else if ( i%2 == 0 && j%2 == 1 && k%2 == 0 ) state = 2; 		
				else if ( i%2 == 1 && j%2 == 1 && k%2 == 0 ) state = 3;
				else if ( i%2 == 0 && j%2 == 0 && k%2 == 1 ) state = 4;
				else if ( i%2 == 1 && j%2 == 0 && k%2 == 1 ) state = 5; 		
				else if ( i%2 == 0 && j%2 == 1 && k%2 == 1 ) state = 6; 		
				else if ( i%2 == 1 && j%2 == 1 && k%2 == 1 ) state = 7;
				switch(state){
				case (0): // These Points are direct transforms from Coarse to Fine
					fineGrid[i][j][k] = coarseGrid[i/2][j/2][k/2];
					break;
				case (1): 
					// Odd Row, Even Columns
					// Average Points Up and DownA
					/***Left off right here with updating the code... ***/
					fineGrid[i][j][k] = (coarseGrid[i/2][j/2][k/2]+coarseGrid[i/2+1][j/2][k/2])/2.0;
					break;
				case (2):
					// Even Row, Odd Columns
					// Same Idea as Case 1
					fineGrid[i][j][k] = (coarseGrid[i/2][j/2][k/2]+coarseGrid[i/2][j/2+1][k/2])/2.0;
					break;
				case (3):
					// Odd Row Odd Columns
					// Average 4 Coarse Points around this point...
					fineGrid[i][j][k] = (coarseGrid[i/2][j/2][k/2]+
							coarseGrid[i/2+1][j/2][k/2]
							+coarseGrid[i/2][j/2+1][k/2]+coarseGrid[i/2+1][j/2+1][k/2])/4.0;
					break;
				case(4):
					fineGrid[i][j][k] = (coarseGrid[i/2][j/2][k/2]+coarseGrid[i/2][j/2][k/2+1])/2.0;
					break;
				case(5):
					fineGrid[i][j][k] = (coarseGrid[i/2][j/2][k/2]+
							coarseGrid[i/2+1][j/2][k/2]
							+coarseGrid[i/2][j/2][k/2+1]+coarseGrid[i/2+1][j/2][k/2+1])/4.0;
					break;
				case(6):
					fineGrid[i][j][k] = (coarseGrid[i/2][j/2][k/2]+
							coarseGrid[i/2][j/2+1][k/2]
							+coarseGrid[i/2][j/2][k/2+1]+coarseGrid[i/2][j/2+1][k/2+1])/4.0;
					break;
				case(7):
					fineGrid[i][j][k] = (coarseGrid[i/2][j/2][k/2]+
							coarseGrid[i/2][j/2][k/2+1]
							+coarseGrid[i/2][j/2+1][k/2]+coarseGrid[i/2][j/2+1][k/2+1]
							+coarseGrid[i/2+1][j/2][k/2]
							+coarseGrid[i/2+1][j/2][k/2+1]
							+coarseGrid[i/2+1][j/2+1][k/2]
							+coarseGrid[i/2+1][j/2+1][k/2+1]
							)/8.0;
					break;

				}
			} 		
		}
	}
}

void restrictToCoarse(int fromDim,double ***dest, double ***origin){
	int i,j,k;
	int N = fromDim/2 + 1;
	// options for Resctrion 
	// 1. Injection
	// 2. Full Weighing
	// 3. Half Weighing
	// USING HERE: half Weighing | Full weighing is way too complicated in 3D

	for(i=1;i<N-1;i++){
		for(j=1;j<N-1;j++){
			for(k=1;k<N-1;k++){
				dest[i][j][k] = 1.0/12 * (
					origin[2*i-1][2*j][2*k] + origin[2*i+1][2*j][2*k] + origin[2*i][2*j-1][k]
					+ origin[2*i][2*j+1][k] + origin[2*i][2*j][2*k-1] + origin[2*i][2*j][2*k+1]
					+ origin[2*i][2*j][2*k]*6.0
				); 

			}
		}
	}
}


void computeResidual(int dim, double ***u, double ***res, double ***rho){
	int i,j,k, N=dim;
	double h = 1.0/(N + 1.0);
	for(i=1;i<dim-1;i++){
	        for(j=1; j<dim-1; j++){
			for(k=1;k<dim-1;k++){
				res[i][j][k] = (rho[i][j][k]
					 - (6.0*u[i][j][k]
						-u[i+1][j][k]-u[i-1][j][k]
						-u[i][j+1][k]-u[i][j-1][k]
						-u[i][j][k+1]-u[i][j][k-1])/h/h);
			}
		}
	}
}

double residualNorm(int dim, double ***u, double ***res, double ***rho){
	int i,j,k, N=dim;
	double h = 1.0/(N + 1.0);

	double resD=0, temp=0;

	for(i=1;i<dim-1;i++){
	        for(j=1; j<dim-1; j++){
			for(k=1;k<dim-1;k++){
		//		printf("Calculating Residual Norm\n (i,j,k) = (%d,%d,%d)\n",i,j,k);
				temp = (rho[i][j][k]-(6.0*u[i][j][k]-u[i+1][j][k]-u[i-1][j][k]-u[i][j+1][k]-u[i][j-1][k]-u[i][j][k+1]-u[i][j][k-1])/h/h); 
				resD += temp * temp;
		//		printf("Current Resd=%f\n",resD);
			}
		}
	}
	return sqrt(resD);

}

// Do GS relaxation on u with RHS = **rhs....
// Not using Red/Black sweeps, but maybe implement in future 
// For better convergence...
void relax(int iterations, double ***u,double ***rhs, int dim){
	relaxCount++;
	int N = dim;
	double h = 1.0/(N+1.0);
	double w = 6.0 / 7.0;
	int i,j,k, iter;
	//double **unew = myMalloc(N);

	// Doing GS here
	// For Parallel MG -> Consider doing Jacobi Instead
	for (iter=0; iter<iterations; iter++){
		for (i=1; i<dim-1; i++){
			for (j=1; j < dim-1; j++){
				for(k=1;k<dim-1;k++){
	//				if (i == N/2 && j == N/2 && k == N/2 )continue;
					u[i][j][k] = u[i][j][k] + w * (
							   u[i-1][j][k] + u[i+1][j][k]
							 + u[i][j-1][k] + u[i][j+1][k]
							 + u[i][j][k+1] + u[i][j][k-1]
							 + rhs[i][j][k]*h*h - 6.0*u[i][j][k]
							)/6.0;
				}
			}
		}
	}
}
