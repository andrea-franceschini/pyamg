#include "inl_blas1.h"
#include "cblas.h" // to use: dgemv

double fsai_ScalFact(const int mrow, const double a22, const double *const a21,
                     const double *const a12, const double *const a11, 
                     const double *const fl, const double *const fu, double *vscr){

   double fac = a22;
   if (mrow > 0){
      fac -= inl_ddot(mrow,a12,1,fl,1);
      fac -= inl_ddot(mrow,a21,1,fu,1);
      cblas_dgemv(CblasRowMajor,CblasNoTrans,mrow,mrow,1.0,a11,mrow,fu,1,0.0,vscr,1);
      fac += inl_ddot(mrow,vscr,1,fl,1);
   }
   return fac;

}
