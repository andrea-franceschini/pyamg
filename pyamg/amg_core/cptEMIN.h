#include <cstdio>
#include <omp.h>

#include "EMIN_ImpProl.h"

// MAIN PROGRAM
template <class I, class R>
int cptEMIN( I itmax,
             R tol,
             R condmax,
             I precType,
             I nn,
             I nc,
             I ntv,
             I const fcnodes[], int const fcnodes_size,
             I const iat_A[], int const iat_A_size,
             I const ja_A[], int const ja_A_size,
             R const coef_A[], int const coef_A_size,
             I const iat_P0[], int const iat_P0_size,
             I const ja_P0[], int const ja_P0_size,
             R const coef_P0[], int const coef_P0_size,
             R const TV_buff[], int const TV_buff_size,
             I const iat_patt[], int const iat_patt_size,
             I const ja_patt[], int const ja_patt_size,
             I iat_Pout[], int const iat_Pout_size,
             I ja_Pout[], int const ja_Pout_size,
             R coef_Pout[], int const coef_Pout_size,
             R info[], int const info_size ) {

   // Init error code
   int ierr = 0;

   // Retrieve number of openMP threads
   int nthreads=1;
   const char* env_p = std::getenv("OMP_NUM_THREADS");
   if(env_p){
      nthreads = atoi(env_p);
   }
   printf("nthreads %d\n",nthreads);
   printf("precType %d\n",precType);

   // Allocate and set TV
   const double **TV = (const double**) malloc(nn*sizeof(double*));
   if (TV == nullptr) return ierr = 1;
   const double *tmp = TV_buff;
   for (int i = 0; i < nn; i++){
      TV[i] = tmp;
      tmp += ntv;
   }

   // --- Minimize Prolongation Energy ---------------------------------------------------
   int nt_A = iat_A[nn];
   int nt_P0 = iat_P0[nn];
   int nt_patt = iat_patt[nn];
   ierr = EMIN_ImpProl(nthreads,itmax,tol,condmax,precType,nn,nc,ntv,nt_A,nt_P0,nt_patt,
                       fcnodes,iat_A,ja_A,coef_A,iat_P0,ja_P0,coef_P0,iat_patt,ja_patt,
                       TV,iat_Pout,ja_Pout,coef_Pout,info);

   return ierr;

}
