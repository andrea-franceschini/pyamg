#include <iostream>
#include <omp.h>
#include <cmath>
#include <chrono>

#if defined PRINT
#define dump true
#else
#define dump false
#endif

#include "EMIN_parm.h"
#include "EMIN_matfree.h"

/*****************************************************************************************
 *
 * Function that starting from an input prolongation and a given non-zero pattern,
 * improves the prolongation through an energy minimization approach while maintaining the
 * near kernel included in the prolongation range.
 *
 * Input:
 *
 * verb:                            verbosity level.
 * np:                              number of openMP threads.
 * itmax:                           number of EnerMinCG iterations.
 * en_tol:                          exit tolerance controlling energy decrease.
 * condmax:                         max conditioning allowed for a B block.
 * prec_type:                       preconditioner for energy minimization
 * nn:                              size of the system matrix.
 * nn_C:                            number of coarse nodes.
 * ntv:                             number of test vectors.
 * nt_A:                            number of non-zeroes of the system matrix.
 * nt_P:                            number of non-zeroes of the input prolongation.
 * nt_patt:                         number of non-zeroes of the non-zero pattern.
 * fcnode:                          F/C indicator.
 * iat_A, ja_A, coef_A:             system matrix.
 * iat_Pin, ja_Pin, coef_Pin:       input prolongation.
 * iat_patt, ja_patt:               prescribed prolongation pattern.
 * TV:                              nn x ntv test space.
 *
 * Output:
 *
 * iat_Pout, ja_Pout, coef_Pout:    output prolongation
 * info:                            array with timings and information on EMIN process
 *                                  info[0]  --> time_prec_K
 *                                  info[1]  --> time_gath_B
 *                                  info[2]  --> time_PCG
 *                                  info[3]  --> time_overhead
 *                                  info[4]  --> time_glob
 *                                  info[5]  --> iter
 *                                  info[6]  --> nnz_J
 *                                  info[7]  --> nnz_Q
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 * 2 ---> error in preconditioner computation
 * 3 ---> error in gather_B
 * 4 ---> error in PCG
 * 5 ---> allocation for final result
 *
 * NOTE: It is assumed that the input pattern does not contain row entries relative
 * to coarse nodes.
 *
 *****************************************************************************************/

int EMIN_ImpProl(const int verb, const int np, const int itmax, const double en_tol,
                 const double condmax, const int prec_type, const int nn, const int nn_C,
                 const int ntv, const int nt_A, const int nt_P, const int nt_patt,
                 const int *fcnode, const int *iat_A, const int *ja_A, const double *coef_A,
                 const int *iat_Pin, const int *ja_Pin, const double *coef_Pin,
                 const int *iat_patt, const int *ja_patt, const double *const *TV,
                 int *iat_Pout, int *ja_Pout, double *coef_Pout, double *info){

   // Init error code
   int ierr = 0;

   // --- Local variables for timing -----------------------------------------------------
   std::chrono::time_point<std::chrono::system_clock> glob_start, glob_end;
   std::chrono::duration<double> elaps_sec;

   //---GLOBAL START-----------------------------
   glob_start = std::chrono::system_clock::now();
   //--------------------------------------------

   ierr = EMIN_matfree(verb,np,itmax,en_tol,condmax,prec_type,nn,nn_C,ntv,nt_A,nt_P,
                       nt_patt,fcnode,iat_A,ja_A,coef_A,iat_Pin,ja_Pin,coef_Pin,
                       iat_patt,ja_patt,TV,iat_Pout,ja_Pout,coef_Pout,info);

   //---GLOBAL END---------------------------------------------------------------
   glob_end = std::chrono::system_clock::now();
   elaps_sec = glob_end - glob_start;
   double time_glob =  elaps_sec.count();
   double time_overhead = time_glob - info[0] - info[1] - info[2];
   //----------------------------------------------------------------------------

   // Store info
   info[4]  = time_overhead;
   info[5]  = time_glob;

   return ierr;

}
