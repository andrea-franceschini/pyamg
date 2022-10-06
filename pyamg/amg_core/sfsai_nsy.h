#pragma once
//#include <iostream>  // to use: cout,endl
//#include <stdlib.h>  // to use: exit
//#include <fstream>   // to use: ifstream,ofstream
//#include <math.h>
//#include <libgen.h>  // to use: basename
//#include <chrono>    // to use: system_clock,duration

#include "mk_pattern.h"
#include "cpt_nsy_sfsai_coef.h"
#include "compress_nsy_sfsai.h"
#include "transp.h"

// MAIN PROGRAM
template <class I, class R>
int sfsai_nsy( I verbosity,
               I kpow,
               I nnzr_max,
               R tau_pref,
               R tau_post,
               I nn_A,
               I nt_A,
               I const iat_A[], int const iat_A_size,
               I const ja_A[], int const ja_A_size,
               R const coef_A[], int const coef_A_size,
               I iat_FL_out[], int const iat_FL_out_size,
               I ja_FL_out[], int const ja_FL_out_size,
               R coef_FL_out[], int const coef_FL_out_size,
               I iat_FU_out[], int const iat_FU_out_size,
               I ja_FU_out[], int const ja_FU_out_size,
               R coef_FU_out[], int const coef_FU_out_size ){

   // Init error code
   int ierr = 0;

   // Retrieve the number of openMP threads
   int nthreads = 1;
   const char* env_p = std::getenv("OMP_NUM_THREADS");
   if(env_p){
      nthreads = atoi(env_p);
   }

   // --- Compute the preconditioner pattern ---------------------------------------------

   int mmax;
   int *iat_FL = nullptr;
   int *ja_FL = nullptr;
   ierr = mk_pattern(verbosity,nthreads,kpow,tau_pref,nnzr_max,nn_A,nt_A,iat_A,ja_A,
                     coef_A,mmax,iat_FL,ja_FL);
   if (ierr != 0) return ierr = 1;

   if (verbosity >= 1){
      printf("Max number of Non-Zeroes per row:                   %3d\n",mmax);
   }

   // --- Compute the preconditioner coefficients ----------------------------------------

   int nt_F_max = iat_FL[nn_A];
   int *nt_FL = (int*) malloc((nthreads+1) * sizeof(int));
   int *nt_FU = (int*) malloc((nthreads+1) * sizeof(int));
   double *coef_FL = (double*) malloc(nt_F_max * sizeof(double));
   double *coef_FU = (double*) malloc(nt_F_max * sizeof(double));
   if (nt_FL == nullptr || nt_FU == nullptr || coef_FL == nullptr ||
       coef_FU == nullptr) return ierr = 2;
   ierr = cpt_nsy_sfsai_coef(nthreads,tau_post,mmax,nn_A,iat_A,ja_A,coef_A,iat_FL,ja_FL,
                             &(nt_FL[1]),&(nt_FU[1]),coef_FL,coef_FU);
   if (ierr != 0) return ierr = 2;

   // Create pointers for FL and FU
   nt_FL[0] = 0;
   nt_FU[0] = 0;
   for (int ip = 0; ip < nthreads; ip++){
      nt_FL[ip+1] += nt_FL[ip];
      nt_FU[ip+1] += nt_FU[ip];
   }
   if (verbosity >= 1){
      double nn_A_d = static_cast<double>(nn_A);
      double avg_nnzr_FL = static_cast<double>(nt_FL[nthreads]) / nn_A_d;
      double avg_nnzr_FU = static_cast<double>(nt_FU[nthreads]) / nn_A_d;
      printf("Final avg nnzr for FL:                       %10.2f\n",avg_nnzr_FL);
      printf("Final avg nnzr for FU:                       %10.2f\n",avg_nnzr_FU);
   }

   // --- Compress the preconditioner ----------------------------------------------------

   int *iat_FU;
   int *ja_FU;
   ierr = compress_nsy_sfsai(nthreads,nn_A,nt_FL,nt_FU,iat_FL,ja_FL,iat_FU,ja_FU,coef_FL,
                             coef_FU);
   if (ierr != 0) return ierr = 3;
   free(nt_FL);
   free(nt_FU);

   // --- Transpose FU -------------------------------------------------------------------

   int *iat_tmp;
   int *ja_tmp;
   double *coef_tmp;
   ierr = transp_csrmat(nthreads,nn_A,nn_A,iat_FU,ja_FU,coef_FU,iat_tmp,ja_tmp,coef_tmp);
   if (ierr != 0) return ierr = 4;

   // Swap pointers
   free(iat_FU);
   free(ja_FU);
   free(coef_FU);
   iat_FU = iat_tmp;
   ja_FU = ja_tmp;
   coef_FU = coef_tmp;

   // Copy the preconditioner in the output
   #pragma omp parallel for num_threads(nthreads)
   for (int i = 0; i <= nn_A; i++){
      iat_FL_out[i] = iat_FL[i];
      iat_FU_out[i] = iat_FU[i];
   }
   #pragma omp parallel for num_threads(nthreads)
   for (int i = 0; i < iat_FL[nn_A]; i++){
      ja_FL_out[i] = ja_FL[i];
      coef_FL_out[i] = coef_FL[i];
   }
   #pragma omp parallel for num_threads(nthreads)
   for (int i = 0; i < iat_FU[nn_A]; i++){
      ja_FU_out[i] = ja_FU[i];
      coef_FU_out[i] = coef_FU[i];
   }

   // --- Delete scratch vectors ---------------------------------------------------------

   free(iat_FL);
   free(ja_FL);
   free(coef_FL);
   free(iat_FU);
   free(ja_FU);
   free(coef_FU);

   return ierr;

}
