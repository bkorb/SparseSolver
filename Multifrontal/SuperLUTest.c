/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file superlu.c
 * \brief a small 5x5 example
 * 
 * <pre>
 * * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 * </pre>
 */
#include "slu_ddefs.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

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

int main(int argc, char *argv[])
{
/*
 * Purpose
 * =======
 * 
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * Users' Guide to illustrate how to call a SuperLU routine, and the
 * matrix data structures used by SuperLU.
 *
 */
    for(int N = 10; N<=200; N+=10){
        clock_t t = 0;
        int SMax = 10;
        for(int S=0; S<SMax; S++){
            SuperMatrix A, L, U, B;
            //double   *a, *rhs;
            //double   s, u, p, e, r, l;
            //int      *asub, *xa;
            int      *perm_r; /* row permutations from partial pivoting */
            int      *perm_c; /* column permutation vector */
            int      nrhs, info, i, permc_spec;
            superlu_options_t options;
            SuperLUStat_t stat;

            /* Initialize matrix A. */
            /*m = n = 5;
            nnz = 12;
            if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
            if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
            if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
            s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
            a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
            a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
            asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
            asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
            asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
            xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;*/

            int nz = 5*(N-2)*(N-2) + 16*(N-2) + 12;
            int *xa = (int*) malloc((N*N+1) * sizeof(int));
            int *asub = (int*) malloc(nz * sizeof(int));
            double *a = (double*) malloc(nz * sizeof(double));
            double *b = (double*) malloc(N*N*sizeof(double));
            nz = setup_matrix(N, xa, asub, a);
            setup_vector(N, b);
            double x[N*N];

            /* Create matrix A in the format expected by SuperLU. */
            dCreate_CompCol_Matrix(&A, N*N, N*N, nz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
            
            /* Create right-hand side matrix B. */
            //nrhs = 1;
            //if ( !(rhs = doubleMalloc(N*N * nrhs)) ) ABORT("Malloc fails for rhs[].");
            //for (i = 0; i < N*N; ++i) rhs[i] = 1.0;
            dCreate_Dense_Matrix(&B, N*N, 1, b, N*N, SLU_DN, SLU_D, SLU_GE);

            if ( !(perm_r = intMalloc(N*N)) ) ABORT("Malloc fails for perm_r[].");
            if ( !(perm_c = intMalloc(N*N)) ) ABORT("Malloc fails for perm_c[].");

            /* Set the default input options. */
            set_default_options(&options);
            //options.ColPerm = NATURAL;
            options.SymmetricMode = YES;

            /* Initialize the statistics variables. */
            StatInit(&stat);

            clock_t before = clock();

            /* Solve the linear system. */
            dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

            clock_t diff = clock() - before;
            t += diff;
            
            //dPrint_CompCol_Matrix("A", &A);
            //dPrint_CompCol_Matrix("U", &U);
            //dPrint_SuperNode_Matrix("L", &L);
            //print_int_vec("\nperm_r", N*N, perm_r);

            /* De-allocate storage */
            SUPERLU_FREE (b);
            SUPERLU_FREE (perm_r);
            SUPERLU_FREE (perm_c);
            Destroy_CompCol_Matrix(&A);
            Destroy_SuperMatrix_Store(&B);
            Destroy_SuperNode_Matrix(&L);
            Destroy_CompCol_Matrix(&U);
            StatFree(&stat);
        }
        //printf("Size: %d x %d \n", N, N);
        printf("%lu \n", 1000000*t / CLOCKS_PER_SEC / SMax);
    }
}
