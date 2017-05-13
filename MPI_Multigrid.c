/*
 *	Joe Zuhusky
 *	MPI Implementation of 2D Multigrid-Method
 *	for 2D Poission Equation

 */

#include "util.h"
#include "helper.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

#define v1 3
#define v2 3

int relaxCount=0;
int pGlobal;

void loadBuffers(double *topOut, double *botOut, double *leftOut, double *rightOut, double **u, int N, int p);

void interpolateToFine(int fromDim, double **coarseGrid, double **fineGrid, int rank, int p,
	double *topGhost, double *botGhost, double *leftGhost, double *rightGhost);

void restrictToCoarse(int fromDim,double **dest, double **origin, int rank, int p,
	double *topGhost, double *botGhost, double *leftGhost, double *rightGhost);

void computeResidual(int dim, double **u, double **res, double **rho,
	double *topGhost, double *botGhost, double *leftGhost, double *rightGhost
	);

double residualNorm(int dim, double **u, double **res, double **rho);

void relax(int iterations, double **u,double **rhs, int dim, int rank, int p,
		 double *topGhost, double *botGhost, double *leftGhost, double *rightGhost);

void doSwaps(double *topOut, double *botOut, double *leftOut, double *rightOut,
		double *topGhost, double *botGhost, double *leftGhost, double *rightGhost, int N);

int main(int argc, char **argv){

	/*MPI Schtuff */
	MPI_Init(&argc,&argv);
	int rank, p, coarsestMesh, procsPerRow;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	procsPerRow = sqrt(p);
	pGlobal = p;

	/* Assert Params*/
	if ( argc < 2 ){
		if ( rank == 0 ){
		printf("Please enter the correct number of parameters\n");
		printf("[Dimension of Problem (power of 2)] | [Number Of Grids]\n");
		}
		exit(1);
	}
	/*Check to make sure parallel Job requested */
	if ( p == 1 ){
		printf("This MPI Job requires more than one MPI Task\n");
		MPI_Abort(MPI_COMM_WORLD,-1);
	}

	// Data Arrays for Pointers to Grids & more vars
	int numVcycles, i, j, k, N, numberOfGrids, n2,C;
	N             = atoi(argv[1]); // Input Initial Grid Size (Use power of two)
	N = N / procsPerRow;
	numberOfGrids = floor(log2(N));
	if ( rank == 0 )
		printf("Working with %d Points per Dim, %d points Per Dim per Proc\n",atoi(argv[1]),N);
	if ( argc > 2 ){
		int g = atoi(argv[2]);
		if ( g > numberOfGrids || g == 0 ){
			if ( rank == 0 ){
				printf("Please a valid number of Grids, should be <= LogN\n");
				exit(0);
			}
		}else{
			numberOfGrids = g;
			if ( rank == 0 ){
				printf("Manual Set # Grids = %d\n",numberOfGrids);
			}
		}
	}


	double **temp, **uPtrs[numberOfGrids], **defects[numberOfGrids], **rhs[numberOfGrids], **tempDefect;
	double **current; // A temp pointer
	double r_norm, initNorm, gNorm, conv_tol = 0.00000001;

	// Declare Buffers for 2D Communication
	double *topGhosts[numberOfGrids], *botGhosts[numberOfGrids];
	double *rightGhosts[numberOfGrids], *leftGhosts[numberOfGrids];
	double *rightOut[numberOfGrids], *leftOut[numberOfGrids];
	double *botOut[numberOfGrids], *topOut[numberOfGrids];

	timestamp_type time1,time2;

	if ( powerOfTwo(N) == 0 ){
		printf("Please use a Dim that is a power of 2\n");
		exit(0);
	}else if ( pow(2,numberOfGrids) > N ){
		printf("Please use a smaller number of grids, #grids <= log2(N)\n");
		exit(0);
	}

	// Increase Dim by 1 for multigrid
	N   = N + 1;
	n2 = N;

	// Preprocessing -> Allocate Space for all of the Grids/Defects/Buffers
	for(i=0;i<numberOfGrids;i++){
		uPtrs[i]   = myMalloc(n2);
		rhs[i]     = myMalloc(n2);
		defects[i] = myMalloc(n2);
		// Officially Allocate Space for Buffers
		topGhosts[i] = (double*)calloc(n2,sizeof(double));
		botGhosts[i] = (double*)calloc(n2,sizeof(double));
		leftGhosts[i] = (double*)calloc(n2,sizeof(double));
		rightGhosts[i] = (double*)calloc(n2,sizeof(double));
		topOut[i] = (double*)calloc(n2,sizeof(double));
		botOut[i] = (double*)calloc(n2,sizeof(double));
		leftOut[i] = (double*)calloc(n2,sizeof(double));
		rightOut[i] = (double*)calloc(n2,sizeof(double));

		zeroOut(n2,uPtrs[i]);
		zeroOut(n2,defects[i]);
		zeroOut(n2,rhs[i]);
		n2 = n2 / 2 + 1;
	}

	// Initialize Source Function - Simple one For now
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			rhs[0][i][j] = 1.0;
		}
	}

	// For each C, do a V-Cycle
	get_timestamp(&time1);	
	// Initialize Norm
	k=0;
	computeResidual(N,uPtrs[k],defects[k],rhs[k],topGhosts[k],botGhosts[k],leftGhosts[k],rightGhosts[k]);
	r_norm   = residualNorm(N,uPtrs[0],defects[0], rhs[0]);
	// Gather all norms
	MPI_Reduce(&r_norm,&gNorm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	gNorm    = sqrt(gNorm);
	MPI_Bcast(&gNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	initNorm = gNorm;

	while ( gNorm / initNorm > conv_tol) {
		if (rank == 0 ){
			printf("Norm=%0.8f -> %0.8f %% of init\n",gNorm,gNorm/initNorm * 100 );	
		}
		// Down Stroke of V Cycle
		for (k=0; k<numberOfGrids-1;k++){

			relax(v1,uPtrs[k],rhs[k], N, rank, p,topGhosts[k],botGhosts[k],leftGhosts[k],rightGhosts[k]);
			computeResidual(N,uPtrs[k],defects[k],rhs[k],topGhosts[k],botGhosts[k],leftGhosts[k],rightGhosts[k]);
			restrictToCoarse(N,rhs[k+1],defects[k], rank, p,topGhosts[k],botGhosts[k],leftGhosts[k],rightGhosts[k]);
			N = N / 2 + 1;
			zeroOut(N,uPtrs[k+1]);
		}

		coarsestMesh = N + procsPerRow*(N-1);
		// Exact Solve
		relax(50*v1,uPtrs[numberOfGrids-1],rhs[numberOfGrids-1], N, 
			rank, p,topGhosts[numberOfGrids-1],botGhosts[numberOfGrids-1],leftGhosts[numberOfGrids-1],rightGhosts[numberOfGrids-1]);

		// Up Stroke of V-Cycle
		for (k=numberOfGrids-1; k > 0;k--){

			current = myMalloc(2*N-1); // Some Temp Memory to holder interpolated value...
                        interpolateToFine(N,uPtrs[k],current, rank, p, topGhosts[k],botGhosts[k],leftGhosts[k],rightGhosts[k]);
			N = 2*N-1;
			addMatrices(uPtrs[k-1],current,N);

			relax(v2,uPtrs[k-1],rhs[k-1], N, rank, p,topGhosts[k-1],botGhosts[k-1],leftGhosts[k-1],rightGhosts[k-1]);
			myFree(current,N); // Find a way to get rid of this...

		}
		k=0;
		computeResidual(N,uPtrs[k],defects[k],rhs[k],topGhosts[k],botGhosts[k],leftGhosts[k],rightGhosts[k]);
		r_norm = residualNorm(N,uPtrs[0],defects[0], rhs[0]);
		MPI_Reduce(&r_norm,&gNorm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		gNorm  = sqrt(gNorm);
		MPI_Bcast(&gNorm,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//	if (rank == 0 ) printf("***************\n");
	} // END V-Cycle Loop
	get_timestamp(&time2);

	//printMatPlot(N,uPtrs[0],rank,sqrt(p));
	
	/*Clean Up and Finish */
	// Preprocessing -> Allocate Space for all of the Grids/Defects/Buffers
	for(i=0;i<numberOfGrids;i++){
		free(uPtrs[i]);
		free(rhs[i]);
		free(defects[i] );
		// Officially Allocate Space for Buffers
		free(topGhosts[i]);
		free(botGhosts[i] );
		free(leftGhosts[i]);
		free(rightGhosts[i]);
		free(topOut[i]);
		free(botOut[i]);
		free(leftOut[i]);
		free(rightOut[i]);
	}

	// Print Final MPI & Timing Results
	double elapsed = timestamp_diff_in_seconds(time1,time2), elapsedSum=0.0;
	MPI_Reduce(&elapsed,&elapsedSum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	elapsedSum = elapsedSum / p;
	if (rank == 0 ){
		printf("\n\n============================================\n\n");
		printf("Avg Time for Computation per Processor is %f seconds.\n", elapsedSum);
		printf("Number of Grids: %d\n",numberOfGrids);
		printf("Coarsest Mesh: N=%d\n",coarsestMesh);
		printf("\n============================================\n\n");
        }
	MPI_Finalize();
	return 0;


}

void loadBuffers(double *topOut, double *botOut, double *leftOut, double *rightOut, double **u, int N, int p ){

	int i, procsPerRow = sqrt(p);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if ( rank%procsPerRow > 0  && rank >= procsPerRow){
		for(i=0;i<N;i++){
			leftOut[i] = u[i][1];
			botOut[i] = u[1][i];
			topOut[i] = u[N-1][i];
			rightOut[i] = u[i][N-1];
		}
	}else if (rank % procsPerRow == 0){
		for(i=0;i<N;i++){
			leftOut[i] = u[i][0];
			botOut[i] = u[1][i];
			topOut[i] = u[N-1][i];
			rightOut[i] = u[i][N-1];
		}
	}else if (rank < procsPerRow ){
		for(i=0;i<N;i++){
			leftOut[i] = u[i][1];
			botOut[i] = u[0][i];
			topOut[i] = u[N-1][i];
			rightOut[i] = u[i][N-1];
		}
	}else{
		printf("ERR -> SOME PROCESSOR OUTBUFFERS NOT LOADED\n"); // Blow up...
		exit(0);
	}
}

void doSwaps(double *topOut, double *botOut, double *leftOut, double *rightOut,
		double *topGhost, double *botGhost, double *leftGhost, double *rightGhost, int N){

	int rank, p, procsPerRow;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	MPI_Status status1,status2;

	/* For 2D case, processors per row is the sqrt of the total number of procs
 	 * Assuming that the Number of Procs is a power of 4
 	 * As well, we are partitioning the Mesh into smaller squares
 	 * In the future perhaps we could look at different Partitioning methods
 	 * */
	procsPerRow = sqrt(p);
	/**NOW SWAPPING TOP AND BOTTOM GHOST VECTORS **/
	// If processor is on Bottom Row -> Only Send/Recv Top buffer
	if ( rank < procsPerRow ){
		MPI_Send(topOut,N,MPI_DOUBLE,rank+procsPerRow,0, MPI_COMM_WORLD);
		MPI_Recv(topGhost,N,MPI_DOUBLE,rank+procsPerRow,0,MPI_COMM_WORLD, &status1);
	}
	// Processor On Top Row -> Only Recv/Send bottom
	else if ( (rank+procsPerRow) > (p-1) ){
		MPI_Recv(botGhost,N, MPI_DOUBLE,rank-procsPerRow,0, MPI_COMM_WORLD, &status2);   // Get Bottom
		MPI_Send(botOut,  N, MPI_DOUBLE,rank-procsPerRow,0, MPI_COMM_WORLD);             // Send Bottom
	}
	// Somewhere in the Middle -> Get/send Bottom Ghost, then Send / Recv Top!!
	else{
		MPI_Recv(botGhost,N, MPI_DOUBLE,rank-procsPerRow,0, MPI_COMM_WORLD, &status2);   // Get Bottom
		MPI_Send(botOut,  N, MPI_DOUBLE,rank-procsPerRow,0, MPI_COMM_WORLD);             // Send Bottom
		MPI_Send(topOut,N,MPI_DOUBLE,rank+procsPerRow,0, MPI_COMM_WORLD);
		MPI_Recv(topGhost,N,MPI_DOUBLE,rank+procsPerRow,0,MPI_COMM_WORLD,&status1);
	} 
	/** NOW SWAPPING LEFT AND RIGHT GHOST VECTORS!!!**/
	// If processor is on Leftmost Column -> Send and Recv Right ghost
	if ( rank % procsPerRow == 0 ){
		MPI_Send(rightOut,N,MPI_DOUBLE,rank+1,0, MPI_COMM_WORLD);
		MPI_Recv(rightGhost,N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,&status1);
	}
	// Processor is on Rightmost Column -> Recv / Send left Ghost
	else if ( (rank+1)%procsPerRow == 0 ){
		MPI_Recv(leftGhost,N, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD, &status2);
		MPI_Send(leftOut,  N, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD);
	}
	else{
		MPI_Recv(leftGhost,N, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD, &status2);
		MPI_Send(leftOut,  N, MPI_DOUBLE,rank-1,0, MPI_COMM_WORLD);
		MPI_Send(rightOut,N,MPI_DOUBLE,rank+1,0, MPI_COMM_WORLD);
		MPI_Recv(rightGhost,N,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,&status1);
	}
}


void interpolateToFine(int fromDim, double **coarseGrid, double **fineGrid, int rank, int p,
	double *topGhost, double *botGhost, double *leftGhost, double *rightGhost){
	int i,j;
	int N = fromDim*2 - 1;
	int state;
	// Using Bi-linear Interpolation
	int procsPerRow = sqrt(p);
	double *topOut = (double*)calloc(N,sizeof(double));
	double *botOut = (double*)calloc(N,sizeof(double));
	double *leftOut = (double*)calloc(N,sizeof(double));
	double *rightOut = (double*)calloc(N,sizeof(double));
	loadBuffers(topOut,botOut,leftOut,rightOut,coarseGrid,fromDim,p);
	doSwaps(topOut,botOut,leftOut,rightOut,topGhost,botGhost,leftGhost,rightGhost,fromDim);
	free(topOut);
	free(botOut);
	free(leftOut);
	free(rightOut);
		
	// Copy bottom and Left Ghosts into Original Array
	for(i=0;i<fromDim;i++){
		coarseGrid[0][i] = botGhost[i];
		coarseGrid[i][0] = leftGhost[i];
	}
	// Doing Column Number 1
	for(i=1;i<N-1;i++){
		if ( i%2==0 ){
			fineGrid[N-1][i] = coarseGrid[(N-1)/2][i/2];
			fineGrid[i][N-1] = coarseGrid[i/2][(N-1)/2];
		}else{
			fineGrid[N-1][i] = (coarseGrid[(N-1)/2][i/2] + coarseGrid[(N-1)/2][i/2+1])/2.0;
			fineGrid[i][N-1] = (coarseGrid[i/2][(N-1)/2] + coarseGrid[i/2+1][(N-1)/2])/2.0;
		}
	}
	fineGrid[N-1][N-1] = (topGhost[(N-1)/2]+rightGhost[(N-1)/2]+coarseGrid[(N-1)/2-1][(N-1)/2]+coarseGrid[(N-1)/2][(N-1)/2-1])/4.0;
	fineGrid[N-1][1] = (coarseGrid[(N-1)/2][1] + leftGhost[(N-1)/2] )/2.0;
	fineGrid[1][N-1] = (coarseGrid[1][(N-1)/2] + botGhost[(N-1)/2] )/2.0;

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
					fineGrid[i][j] = (coarseGrid[i/2][j/2]+coarseGrid[i/2+1][j/2])/2.0;
					break;
				case (2):
					// Even Row, Odd Columns
					fineGrid[i][j] = (coarseGrid[i/2][j/2]+coarseGrid[i/2][j/2+1])/2.0;
					break;
				case (3):// Odd Row Odd Columns, Avg 4 corners
					fineGrid[i][j] = (coarseGrid[i/2][j/2]+coarseGrid[i/2+1][j/2]+coarseGrid[i/2][j/2+1]+coarseGrid[i/2+1][j/2+1])/4.0;
					break;
			} 	
		}
	}
}

void restrictToCoarse(int fromDim,double **dest, double **origin, int rank, int p,
	double *topGhost, double *botGhost, double *leftGhost, double *rightGhost){
	int i,j;
	int procsPerRow = sqrt(p);
	int N = fromDim;
	// USING HERE: Half!! Weighing
	double *topOut = (double*)calloc(N,sizeof(double));
	double *botOut = (double*)calloc(N,sizeof(double));
	double *leftOut = (double*)calloc(N,sizeof(double));
	double *rightOut = (double*)calloc(N,sizeof(double));
	loadBuffers(topOut,botOut,leftOut,rightOut,origin,fromDim,p);
	doSwaps(topOut,botOut,leftOut,rightOut,topGhost,botGhost,leftGhost,rightGhost,fromDim);
	N = fromDim/2 + 1;
	free(topOut);
	free(botOut);
	free(leftOut);
	free(rightOut);

	for(i=1;i<N-1;i++){
		for(j=1;j<N-1;j++){
			dest[i][j] = 1.0/8.0* (
				origin[2*i][2*j+1]+origin[2*i][2*j-1]+origin[2*i+1][2*j]+origin[2*i-1][2*j]+4.0*origin[2*i][2*j]); 
		}
		// Right Edge Coarsen
		dest[i][N-1] = (4.0*origin[i][N-1]+ origin[2*i-1][2*(N-1)] + origin[2*i+1][2*(N-1)] + origin[2*i][2*(N-1)-1] + rightGhost[2*i])/8.0;
		// Top Edge Coarsen
		dest[N-1][i] = (4.0*origin[N-1][i]+ origin[2*(N-1)][2*i-1]+ origin[2*(N-1)][2*i+1]+origin[2*(N-1)-1][2*i]+topGhost[2*i] )/8.0;
	}
	dest[N-1][N-1] = (4.0*origin[N-1][N-1]+ rightGhost[2*i] + topGhost[2*i] + origin[2*(N-1)-1][2*(N-1)] + origin[2*(N-1)][2*(N-1)-1] ) / 8.0;
}

void computeResidual(int dim, double **u, double **res, double **rho,double *topGhost, double *botGhost, double *leftGhost, double *rightGhost){

	int i,j, N=dim;
	int ppr = sqrt(pGlobal);
	double h = 1.0/(N - 1)/ppr;
	double *topOut = (double*)calloc(N,sizeof(double));
	double *botOut = (double*)calloc(N,sizeof(double));
	double *leftOut = (double*)calloc(N,sizeof(double));
	double *rightOut = (double*)calloc(N,sizeof(double));
	loadBuffers(topOut,botOut,leftOut,rightOut,u,N,pGlobal);
	doSwaps(topOut,botOut,leftOut,rightOut,topGhost,botGhost,leftGhost,rightGhost,N);
	free(topOut);
	free(botOut);
	free(leftOut);
	free(rightOut);
	for(i=0;i<N;i++){
		u[0][i] = botGhost[i];
		u[i][0] = leftGhost[i];
	}
	for(i=1;i<dim-1;i++){
	        for(j=1; j<dim-1; j++){
			res[i][j] = (rho[i][j] - (4.0*u[i][j]-u[i+1][j]-u[i-1][j]-u[i][j+1]-u[i][j-1])/h/h);
		}
		res[dim-1][i] = (rho[dim-1][i] - (4.0 * u[dim-1][i]- u[dim-2][i]-u[dim-1][i+1]-u[dim-1][i-1] - topGhost[i])/h/h);
		res[i][N-1] = (rho[i][N-1] - (4.0*u[i][N-1] - rightGhost[i] - u[i+1][N-1]-u[i-1][N-1]-u[i][N-2])/h/h);
	}
	res[N-1][N-1] = (rho[N-1][N-1] - (4.0*u[N-1][N-1] - rightGhost[N-1] - topGhost[N-1]-u[N-2][N-1]-u[N-1][N-2])/h/h);
	
}

double residualNorm(int dim, double **u, double **res, double **rho) {
	int i,j, N=dim;
	double resD=0, temp=0;
	for(i=1;i<dim-1;i++){
	        for(j=1; j<dim-1; j++){
			resD += res[i][j]*res[i][j];
		}
	}
	return resD;

}

void relax(int iterations, double **u,double **rhs, int dim, int rank, int p,
		 double *topGhost, double *botGhost, double *leftGhost, double *rightGhost ){
	int i,j, iter, procsPerRow, N=dim;
	double h, w = 1.0;
	double *topOut = (double*)calloc(N,sizeof(double));
	double *botOut = (double*)calloc(N,sizeof(double));
	double *leftOut = (double*)calloc(N,sizeof(double));
	double *rightOut = (double*)calloc(N,sizeof(double));
	procsPerRow = sqrt(p);
	h = 1.0 / (N - 1.0) / procsPerRow;

	for(iter=0;iter<iterations;iter++){

		loadBuffers(topOut,botOut,leftOut,rightOut,u,N,p);
		doSwaps(topOut,botOut,leftOut,rightOut,topGhost,botGhost,leftGhost,rightGhost,N);
		
		// Copy bottom and Left Ghosts into Original Array
		for(i=0;i<N;i++){
			if ( rank % procsPerRow > 0 )
				u[i][0] = leftGhost[i];
			if ( rank >= procsPerRow )
				u[0][i] = botGhost[i];
		}
		for(i=1;i<N-1;i++){
			for(j=1;j<N-1;j++){
				u[i][j] = u[i][j] + w * (u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1]+rhs[i][j]*h*h-4.0*u[i][j])/4.0;
			}
		}
		// Relax on Edges
		for(i=1;i<N-1;i++){
			if ( rank < (p-procsPerRow) ){
				u[N-1][i]=u[N-1][i]+w*(u[N-2][i]+u[N-1][i+1]+u[N-1][i-1]+topGhost[i]+rhs[N-1][i]*h*h-4.0*u[N-1][i])/4.0;
				if ( rank % procsPerRow > 0 ){
					u[N-1][1] = u[N-1][1] + w * (u[N-1][2]+u[N-1][0]+u[N-2][1]+topGhost[1]+rhs[N-1][1]*h*h - 4.0*u[N-1][1])/4.0;
				}
			}
			if ( (rank+1) % procsPerRow > 0 ){
				u[i][N-1] = u[i][N-1] + w * (u[i-1][N-1]+u[i+1][N-1]+u[i][N-2]+rightGhost[i]+rhs[i][N-1]*h*h - 4.0*u[i][N-1])/4.0;
				if ( rank < (p-procsPerRow) ){
					u[N-1][N-1] = u[N-1][N-1]+w*(u[N-2][N-1]+u[N-1][N-2]+topGhost[N-1]+rightGhost[N-1]+rhs[N-1][N-1]*h*h-4.0*u[N-1][N-1])/4.0;
				}
			}
		}
		loadBuffers(topOut,botOut,leftOut,rightOut,u,N,p);
		doSwaps(topOut,botOut,leftOut,rightOut,topGhost,botGhost,leftGhost,rightGhost,N);

		for(i=0;i<N;i++){
			if ( rank % procsPerRow > 0 )
				u[i][0] = leftGhost[i];
			if ( rank >= procsPerRow )
				u[0][i] = botGhost[i];
		}
		for(i=N-2;i>0;i--){
			for(j=N-2;j>0;j--){
				u[i][j] = u[i][j] + w * (u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1]+rhs[i][j]*h*h-4.0*u[i][j])/4.0;
			}
		}
		// Relax on Edges
		for(i=1;i<N-1;i++){
			if ( rank < (p-procsPerRow) ){
				u[N-1][i]=u[N-1][i]+w*(u[N-2][i]+u[N-1][i+1]+u[N-1][i-1]+topGhost[i]+rhs[N-1][i]*h*h-4.0*u[N-1][i])/4.0;
				if ( rank % procsPerRow > 0 ){
					u[N-1][1] = u[N-1][1] + w * (u[N-1][2]+u[N-1][0]+u[N-2][1]+topGhost[1]+rhs[N-1][1]*h*h - 4.0*u[N-1][1])/4.0;
				}
			}
			if ( (rank+1) % procsPerRow > 0 ){
				u[i][N-1] = u[i][N-1] + w * (u[i-1][N-1]+u[i+1][N-1]+u[i][N-2]+rightGhost[i]+rhs[i][N-1]*h*h - 4.0*u[i][N-1])/4.0;
				if ( rank < (p-procsPerRow) ){
					u[N-1][N-1] = u[N-1][N-1]+w*(u[N-2][N-1]+u[N-1][N-2]+topGhost[N-1]+rightGhost[N-1]+rhs[N-1][N-1]*h*h-4.0*u[N-1][N-1])/4.0;
				}
			}
		}
	}
	free(topOut);
	free(botOut);
	free(leftOut);
	free(rightOut);
	
}
