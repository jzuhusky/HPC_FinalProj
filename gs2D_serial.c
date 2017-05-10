
/*
 *	Joe Zuhusky
 *	Serial Implementation of 2D Multigrid-Method
 *	for 2D Poission Equation
 */

#include "util.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv){

	if ( argc != 4 ){
                printf("Please provide 3 Args: Mesh Size, Max Iterations, Convergence Factor (Residual Decrease by this factor to conv)\n");
                exit(0);
        }

	int print_flag = 1;

	const int N = atoi(argv[1])+2;
	int convergenceFactor = atoi(argv[3]);
	int steps=0;
	int i,j;
	long int max_iter = atoi(argv[2]);

	double **u,h;
	double **rho;

	u = (double**)malloc(sizeof(double*)*N);
	rho = (double**)malloc(sizeof(double*)*N);
	
	for (i=0;i<N;i++){
		*(u+i) = (double*)malloc(sizeof(double)*N);
		*(rho+i)  = (double*)malloc(sizeof(double)*N);
	}

	h = 1.0/(N+1.0); 

	// initialize the grid...
	// and note: boundary conditions are such that u(boundary) == 0.0 
	int f,g;
	for (f=0; f<N; f++){
		for (g=0; g<N; g++){
			u[f][g]=0.0;
			rho[f][g]=1.0;
		}
	}
	// Some Charge Distributions
	//rho[N/2][N/2] = -1.0; // Electron Monopole...

	double temp2, res,initRes, conv_tol=1.0/convergenceFactor;
	/* timing */
  	timestamp_type time1, time2;
  	get_timestamp(&time1);
	int numNotConv = 0;
	int FLOPct=0;
	// Doing GS such that f(x,y) = 1.0 for the entire space...
	while( steps == 0 || conv_tol < res/initRes ){
		res = 0.0;
		for (i=1; i<N-1; i++){
			for (j=1; j < N-1; j++){
				u[i][j] = (u[i-1][j]+u[i+1][j] + u[i][j-1] + u[i][j+1] + rho[i][j]*h*h)/4.0;
				FLOPct += 7;;
			}
		}
		for (i=N-2; i>0; i--){
			for (j=N-2; j>0; j--){
				u[i][j] = (u[i-1][j]+u[i+1][j] + u[i][j-1] + u[i][j+1] + rho[i][j]*h*h)/4.0;
				FLOPct += 7;;
			}
		}
		// Computing Residual...
		for(i=1;i<N-1;i++){
			for(j=1; j<N-1; j++){
				temp2 = (4*u[i][j]-u[i-1][j]-u[i+1][j] - u[i][j-1] - u[i][j+1] - rho[i][j]*h*h);
				res += temp2 * temp2;
				FLOPct += 8;
			}
		}
		res = sqrt(res);
		FLOPct++;
		if ( steps == 0 ){
                        initRes = res;
                }
			printf("%d %f\n",steps, res/initRes);
		steps++;	
		if (steps > max_iter) break;
    	} // end while
  	get_timestamp(&time2);
	double elapsed = timestamp_diff_in_seconds(time1,time2);
	printf("Time elapsed is %f seconds.\n", elapsed);	
	printf("%d Iterations & %0.8f percent of initial residual \n",steps,res/initRes*100.0);
	printf("%f GFlops/s\n",FLOPct/1e9/elapsed);
	
	if (print_flag !=0 ){
		for (i=0;i<N;i++){
			for (j=0;j<N;j++){
		//		printf("%f %f %f\n",h*i,h*j,u[i][j]);
			}
		}
	}

	// Free all that malloced space...
	for (i=0;i<N;i++){
		free(*(u+i));
	}
	free(u);
	return 0;
}
