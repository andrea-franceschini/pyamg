/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 * outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include <stdlib.h>
////////////////////////////////
#include <iostream>
////////////////////////////////

#include "BAMG.h"
#include "BAMG_params.h"
#include "DebEnv.h"

// MAIN PROGRAM
template <class I, class R>
int cptBAMGProl( I const level,
                 I const verbosity,
                 I const itmax_vol,
                 I const dist_min,
                 I const dist_max,
                 I const mmax,
                 R const maxcond,
                 R const maxrownrm,
                 R const tol_vol,
                 R const eps,
                 I const nn_S,
                 I const iat_S[], int const iat_S_size,
                 I const ja_S[], int const ja_S_size,
                 I const ntv,
                 I const fcnodes[], int const fcnodes_size,
                 R const TV[], int const TV_size,
                 I const nn_I,
                 I nt_I,
                 I iat_I[], int const iat_I_size,
                 I ja_I[], int const ja_I_size,
                 R coef_I[], int const coef_I_size,
                 I c_mark[], int const c_mark_size )
{
   // Retrieve number of openMP threads
   int nthreads=1;
   const char* env_p = std::getenv("OMP_NUM_THREADS");
   if(env_p){
      nthreads = atoi(env_p);
   }

   // Init debug
   if (level == 1){
      DebEnv.SetDebEnv(nthreads,"w");
   } else {
      DebEnv.OpenDebugLog("a");
   }
   if (DEBUG && BAMG_DEBUG){
      fprintf(DebEnv.r_logfile,"\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n",level);
      fflush(DebEnv.r_logfile);
      for (int i = 0; i < nthreads; i++){
         fprintf(DebEnv.t_logfile[i],"\n+++++++++++++++ LEVEL %2d +++++++++++++++\n\n",level);
         fflush(DebEnv.t_logfile[i]);
      }
   }

   // Load TV in a 2D buffer
   double **TV_2D = (double**) std::malloc(nn_S*sizeof(double*));
   int kk = 0;
   for (int i = 0; i < nn_S; i++){
      TV_2D[i] = (double*) std::malloc(ntv*sizeof(double));
      for (int j = 0; j < ntv; j++) {
         TV_2D[i][j] = TV[kk];
         kk++;
      }
   }

   int * fcnodesBAMG = (int *) std::malloc(nn_S*sizeof(int));
   int cnt = 0;
   for( int i = 0; i < nn_S; ++i )
   {
      if( fcnodes[i] > 0 )
      {
         fcnodesBAMG[i] = cnt++;
      }
      else
      {
         fcnodesBAMG[i] = -1;
      }
   }

   std::vector< iReg > vec_iat_I;
   std::vector< iReg > vec_ja_I;
   std::vector< rExt > vec_coef_I;
   std::vector< iReg > vec_c_mark;

   // Store input in the structure
   BAMG_params params;
   params.verbosity = verbosity;
   params.itmax_vol = itmax_vol;
   params.dist_min = dist_min;
   params.dist_max = dist_max;
   params.mmax = mmax;
   params.maxcond = maxcond;
   params.maxrownrm = maxrownrm;
   params.tol_vol = tol_vol;
   params.eps = eps;

   // Set nn_L to 0 as it is useless without MPI
   iReg nn_L = 0;
   iReg nn_C = nn_S;
   int ierr = BAMG(params,nthreads,nn_L,nn_C,nn_S,iat_S,ja_S,ntv,fcnodesBAMG,TV_2D,
                   nt_I,vec_iat_I,vec_ja_I,vec_coef_I,vec_c_mark);

   // Close debug log
   DebEnv.CloseDebugLog();

   // Save iat_I, c_mark
   for( int i = 0; i < nn_I; ++i )
   {
      iat_I[i] = vec_iat_I[i];
      c_mark[i] = vec_c_mark[i];
   }
   iat_I[nn_I] = vec_iat_I[nn_I];
   vec_iat_I.resize( 0 );
   vec_c_mark.resize( 0 );
   // Save ja_I, coef_I
   for( int i = 0; i < nt_I+1; ++i )
   {
      ja_I[i] = vec_ja_I[i];
      coef_I[i] = vec_coef_I[i];
   }
   vec_ja_I.resize( 0 );
   vec_coef_I.resize( 0 );

   // Free temporarily allocated arrays
   for (int i = 0; i < nn_S; i++)
   {
      std::free(TV_2D[i]);
   }
   std::free(TV_2D);
   std::free( fcnodesBAMG );

   return ierr;
}
