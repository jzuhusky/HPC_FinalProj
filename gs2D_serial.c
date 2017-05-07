
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

	double temp2, res,initRes;
	/* timing */
  	timestamp_type time1, time2;
  	get_timestamp(&time1);
	int numNotConv = 0;
	// Doing GS such that f(x,y) = 1.0 for the entire space...
	printf("#Gs step number | norm\n");
	while( steps < max_iter ){
		res = 0.0;
		for (i=1; i<N-1; i++){
			for (j=1; j < N-1; j++){
				u[i][j] = (u[i-1][j]+u[i+1][j] + u[i][j-1] + u[i][j+1] + rho[i][j]*h*h)/4.0;
			}
		}
		// Computing Residual...
		for(i=1;i<N-1;i++){
			for(j=1; j<N-1; j++){
				temp2 = (4*u[i][j]-u[i-1][j]-u[i+1][j] - u[i][j-1] - u[i][j+1] - rho[i][j]*h*h);
				res += temp2 * temp2;
			}
		}
		res = sqrt(res);
		printf("%d %f\n",steps, res);
		if ( steps == 0 ){
                        initRes = res;
                }else{
                        if ( initRes / res > convergenceFactor ){
                               break;
                        }
                }
		steps++;	
    	} // end while
  	get_timestamp(&time2);

	printf("%d Iterations & %f percent of initial residual \n",steps,res/initRes*100.0);
	
	if (print_flag !=0 ){
		for (i=0;i<N;i++){
			for (j=0;j<N;j++){
		//		printf("%f %f %f\n",h*i,h*j,u[i][j]);
			}
		}
	}

	double elapsed = timestamp_diff_in_seconds(time1,time2);
//	printf("Time elapsed is %f seconds.\n", elapsed);	
//	printf("%d %d %f\n",N-2,atoi(argv[argc-1]) , elapsed);
	// Free all that malloced space...
	for (i=0;i<N;i++){
		free(*(u+i));
	}
	free(u);
	return 0;
}
