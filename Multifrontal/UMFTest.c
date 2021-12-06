/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) 2005-2012 by Timothy A. Davis,                       */
/* http://www.suitesparse.com. All Rights Reserved.                           */
/* See ../Doc/License.txt for License.                                        */
/* -------------------------------------------------------------------------- */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "umfpack.h"

//int    n = 5 ;
//int    Ap [ ] = {0, 2, 5, 9, 10, 12} ;
//int    Ai [ ] = { 0,  1,  0,   2,  4,  1,  2,  3,   4,  2,  1,  4} ;
//double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
//double b [ ] = {8., 45., -3., 3., 19.} ;
//double x [5] ;

double random_d(){
    return (double)rand() / (double)RAND_MAX ;
}

void setup_vector(int N, double *b){
    for(int i = 0; i<N*N; i++){
	b[i] = random_d();
    }
}

int setup_matrix(int N, int *Ap, int *Ai, double *Ax){
    int nz = 5*(N-2)*(N-2) + 16*(N-2) + 12;
    int index = 0;
    for(int j = 0; j<N; j++){
	for(int i = 0; i<N; i++){
	    Ap[N*j + i] = index;
            if(i>0){
		Ai[index] = N*(i-1) + j;
		Ax[index] = -1;
		index++;
	    }
	    if(j>0){
		Ai[index] = N*i + j-1;
		Ax[index] = -1;
		index++;
	    }
	    Ai[index] = N*i + j;
	    Ax[index] = 4;
	    index++;
	    if(j<N-1){
		Ai[index] = N*i + j+1;
		Ax[index] = -1;
		index++;
	    }
	    if(i<N-1){
		Ai[index] = N*(i+1) + j;
		Ax[index] = -1;
		index++;
	    }
	}	
    }
    Ap[N*N] = index;
    return nz;
}

int main (void)
{
    for(int N = 10; N<=200; N+=10){
    	clock_t t = 0;
    	int SMax = 10;
    	for(int S = 0; S<SMax; S++){
		    int nz = 5*(N-2)*(N-2) + 16*(N-2) + 12;
		    int *Ap = (int*) malloc((N*N+1) * sizeof(int));
		    int *Ai = (int*) malloc(nz * sizeof(int));
		    double *Ax = (double*) malloc(nz * sizeof(double));
			double *b = (double*) malloc(N*N*sizeof(double));
			nz = setup_matrix(N, Ap, Ai, Ax);
			setup_vector(N, b);
			double x[N*N];

			//Start a timer
			clock_t before = clock();

	        double *null = (double *) NULL ;
	        int i ;
	        void *Symbolic, *Numeric ;
	        (void) umfpack_di_symbolic (N*N, N*N, Ap, Ai, Ax, &Symbolic, null, null) ;
	        (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
	        umfpack_di_free_symbolic (&Symbolic) ;
	        (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null) ;
	        umfpack_di_free_numeric (&Numeric) ;

	        clock_t diff = clock() - before;
	        t += diff;
	        //for (i = 0 ; i < N*N ; i++) printf ("x [%d] = %g, ", i, x [i]) ;
	        //printf("\n");
			free(Ap);
			free(Ai);
			free(Ax);
			free(b);
		}
		//printf("Size: %d x %d \n", N, N);
		printf("%lu \n", 1000000*t / CLOCKS_PER_SEC / SMax);
    }
    return (0) ;
}

