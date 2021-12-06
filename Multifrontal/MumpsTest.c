/*
 *
 *  This file is part of MUMPS 5.4.1, released
 *  on Tue Aug  3 09:49:43 UTC 2021
 *
 */
/* Example program using the C interface to the 
 * double real arithmetic version of MUMPS, dmumps_c.
 * We solve the system A x = RHS with
 *   A = diag(1 2) and RHS = [1 4]^T
 * Solution is [1 2]^T */
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"
#include "dmumps_c.h"
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

double random_d(){
    return (double)rand() / (double)RAND_MAX ;
}

void setup_vector(int N, double *b){
    for(int i = 0; i<N*N; i++){
    	b[i] = random_d();
    }
}

int setup_matrix(int N, int *irn, int *jcn, double *a){
    int nz = 5*(N-2)*(N-2) + 16*(N-2) + 12;
    int index = 0;
    for(int i = 0; i<N; i++){
	    for(int j = 0; j<N; j++){
	        if(i>0){
		        irn[index] = N*i+j+1;
		        jcn[index] = N*(i-1)+j+1;
		        a[index] = -1;
		        index++;
	        }
	        if(j>0){
	        	irn[index] = N*i+j+1;
	        	jcn[index] = N*i+j-1+1;
		        a[index] = -1;
		        index++;
	        }
	        irn[index] = N*i+j+1;
	        jcn[index] = N*i+j+1;
	        a[index] = 4;
	        index++;
	        if(j<N-1){
	        	irn[index] = N*i+j+1;
	        	jcn[index] = N*i+j+1+1;
		        a[index] = -1;
		        index++;
	        }
	        if(i<N-1){
	        	irn[index] = N*i+j+1;
	        	jcn[index] = N*(i+1)+j+1;
		        a[index] = -1;
		        index++;
	        }
	    }   
    }
    return nz;
}

int main(int argc, char ** argv){
	int myid, ierr;

	int error = 0;
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	for(int N = 10; N<=200; N+=10){
		clock_t t = 0;
        int SMax = 10;
		for(int S = 0; S<SMax; S++){
		  DMUMPS_STRUC_C id;
		  int n = N*N;
		  int nz = 5*(N-2)*(N-2) + 16*(N-2) + 12;
		  int *irn = (int*) malloc(nz * sizeof(int));
		  int *jcn = (int*) malloc(nz * sizeof(int));
		  double *a = (double*) malloc(nz * sizeof(double));
		  double *b = (double*) malloc(N*N*sizeof(double));
		  setup_matrix(N, irn, jcn, a);
		  //printf("%d\n", N);
		  setup_vector(N, b);
		/* When compiling with -DINTSIZE64, MUMPS_INT is 64-bit but MPI
		   ilp64 versions may still require standard int for C interface. */
		/* MUMPS_INT myid, ierr; */
		  
		  /* Define A and rhs */
		  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
		  id.comm_fortran=USE_COMM_WORLD;
		  id.par=1; id.sym=0;
		  id.job=JOB_INIT;
		  dmumps_c(&id);

		  /* Define the problem on the host */
		  if (myid == 0) {
		    id.n = N*N; id.nnz = nz; id.irn=irn; id.jcn=jcn;
		    id.a = a; id.rhs = b;
		  }

		  /* No outputs */
		  id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0; //id.ICNTL(5)=0; id.ICNTL(18)=0;

		  /* Call the MUMPS package (analyse, factorization and solve). */
		  id.job=6;

		  clock_t before = clock();

		  dmumps_c(&id);

		  clock_t diff = clock() - before;
          t += diff;

		  if (id.infog[0]<0) {
		    printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
		        myid, id.infog[0], id.infog[1]);
		    error = 1;
		  }

		  /* Terminate instance. */
		  id.job=JOB_END;
		  dmumps_c(&id);
		  if (myid == 0) {
		    if (!error) {
		      //printf("Solution is : (%f  %f)\n", b[0],b[1]);
		    } else {
		      //printf("An error has occured, please check error code returned by MUMPS.\n");
		    }
		  }
		}
		//printf("Size: %d x %d \n", N, N);
        printf("%lu \n", 1000000*t / CLOCKS_PER_SEC / SMax);
	}
	ierr = MPI_Finalize();
  return 0;
}
