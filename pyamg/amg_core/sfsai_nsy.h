#pragma once

#include <cstdio>

#include "mk_pattern.h"
#include "cpt_nsy_sfsai_coef.h"
#include "compress_nsy_sfsai.h"

// MAIN PROGRAM
template <class I, class R>
int sfsai_nsy( I kpow,
               I nnzr_max,
               R tau_pref,
               R tau_post,
               I nn_A,
               I nt_A,
               I const iat_A[], int const iat_A_size,
               I const ja_A[], int const ja_A_size,
               R const coef_A[], int const coef_A_size,
               I iat_FL[], int const iat_FL_size,
               I ja_FL[], int const ja_FL_size,
               R coef_FL[], int const coef_FL_size,
               I iat_FU[], int const iat_FU_size,
               I ja_FU[], int const ja_FU_size,
               R coef_FU[], int const coef_FU_size )
{

   int ierr = 0;

   // --- Compute the preconditioner pattern ---------------------------------------------

   int mmax;
   ierr = mk_pattern(kpow,tau_pref,nnzr_max,nn_A,nt_A,iat_A,ja_A,coef_A,mmax,iat_FL,ja_FL_size,ja_FL);
   if (ierr != 0) return ierr;

   // --- Compute the preconditioner coefficients ----------------------------------------

   int nt_FL;
   int nt_FU;
   ierr = cpt_nsy_sfsai_coef(tau_post,mmax,nn_A,iat_A,ja_A,coef_A,iat_FL,ja_FL,
                             nt_FL,nt_FU,coef_FL,coef_FU);
   if (ierr != 0) return ierr;

   // --- Compress the preconditioner (transposing FU) -----------------------------------

   ierr = compress_nsy_sfsai(nn_A,nt_FL,nt_FU,iat_FL,ja_FL,iat_FU,ja_FU,coef_FL,coef_FU);

   return ierr;
}
