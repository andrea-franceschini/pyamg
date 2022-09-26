//#include <iostream>  // to use: cout,endl
#include <cstdio>
#include <cmath>

#include "transp_Patt.h"
#include "merge_row_patt.h"
#include "power_patt.h"

int mk_pattern(const int verb, const int kpow, const double tau_pref, const int nnzr_max,
               const int nn_A, const int nt_A, const int *const iat_A, const int *const ja_A,
               const double *const coef_A, int &mmax, int * iat_Pout, const int ja_Pout_size,
               int * ja_Pout){

   // Init error code
   int ierr = 0;

   // Extract diagonals form A
   double *diag_A = new double [nn_A]();
   if (diag_A == nullptr) return 1;
   for (int i = 0; i < nn_A; i++){
      int ind = iat_A[i];
      int jcol = ja_A[ind];
      //@@@@ HERE WE ASSUME THAT EACH ROW CONTAINS A NON-ZERO DIAGONAL ENTRY
      while (jcol < i){
         ind++;
         jcol = ja_A[ind];
      }
      diag_A[i] = std::abs(coef_A[ind]);
   }

   // Filter A
   int *iat_FA =  new int [nn_A+1]();
   if (iat_FA == nullptr) return 1;
   int *ja_FA  =  new int [nt_A]();
   if (ja_FA == nullptr) return 1;
   int nt_FA = 0;
   for (int i = 0; i < nn_A; i++){
      iat_FA[i] = nt_FA;
      double fac = tau_pref*diag_A[i];
      for (int j = iat_A[i]; j < iat_A[i+1]; j++){
         int jcol = ja_A[j];
         if (pow(std::abs(coef_A[j]),2) > fac*diag_A[jcol]){
            ja_FA[nt_FA] = jcol;
            nt_FA++;
         }
      }
   }
   iat_FA[nn_A] = nt_FA;
   // Delete scratch
   delete [] diag_A;
   /////////////////////////////////////////////////////
   //FILE *ofile = fopen("mat_InitPatt.csr","w"); if (!ofile) exit(1);
   //for (int i = 0; i < nn_A; i++){
   //    for (int j = iat_FA[i]; j < iat_FA[i+1]; j++){
   //       fprintf(ofile,"%10d %10d %1d\n",i+1,ja_FA[j]+1,1);
   //    }
   //}
   //fclose(ofile);
   /////////////////////////////////////////////////////

   double nn_A_d = static_cast<double>(nn_A);
   if (verb >= 2){
      double nt_FA_d = static_cast<double>(nt_FA);
      fprintf(stdout,"Initial Filtered Pattern avg nnzr:           %10.2f\n",nt_FA_d / nn_A_d);
   }

   // Transpose filtered pattern
   int *iat_FT;
   int *ja_FT;
   ierr = transp_Patt(nn_A,nn_A,iat_FA,ja_FA,iat_FT,ja_FT);
   if (ierr != 0) return 2;

   // Merge filtered and filtered transposed patterns
   int nt_F_max = 2*nt_FA;
   int *iat_F = new int [nn_A+1]();
   if (iat_F == nullptr) return 1;
   int *ja_F = new int [nt_F_max]();
   if (ja_F == nullptr) return 1;
   int ind = 0;
   int iend_FA = iat_FA[0];
   int iend_FT = iat_FT[0];
   for (int i = 0; i < nn_A; i++){
      int len_out;

      // Set pointer to output
      iat_F[i] = ind;
      // Set pointer to input 1
      int istrt_FA = iend_FA;
      iend_FA = iat_FA[i+1];
      int len_FA = iend_FA-istrt_FA;
      // Set pointer to input 2
      int istrt_FT = iend_FT;
      iend_FT = iat_FT[i+1];
      int len_FT = iend_FT-istrt_FT;
      // Merge
      merge_row_patt(len_FA,&(ja_FA[istrt_FA]),len_FT,&(ja_FT[istrt_FT]),len_out,&(ja_F[ind]));
      // Update pointer
      ind += len_out;
   }
   iat_F[nn_A] = ind;

   // Delete intermediate scratches
   delete [] iat_FT;
   delete [] ja_FT;
   delete [] iat_FA;
   delete [] ja_FA;

   //int nt_F = ind;
   //double nt_F_d = static_cast<double>(nt_F);
   //fprintf(stdout,"Final Filtered Symmetrized Pattern avg nnzr: %10.2f\n",nt_F_d / nn_A_d);

   /////////////////////////////////////////////////////
   //ofile = fopen("mat_FiltPatt.csr","w"); if (!ofile) exit(1);
   //for (int i = 0; i < nn_A; i++){
   //    for (int j = iat_F[i]; j < iat_F[i+1]; j++){
   //       fprintf(ofile,"%10d %10d %1d\n",i+1,ja_F[j]+1,1);
   //    }
   //}
   //fclose(ofile);
   /////////////////////////////////////////////////////

   // Compute the power of the final filtered pattern
   ierr = power_patt(kpow,nnzr_max,nn_A,iat_F,ja_F,mmax,iat_Pout,ja_Pout_size,ja_Pout);
   if (ierr != 0)
   {
      //if (ierr < 0) std::cout << "Increase pattern size up to " << -ierr << std::endl;
      return 2;
   }

   if (verb >= 2){
      int max_nnzr = 0;
      for (int i = 0; i < nn_A; i++){
         int len = iat_Pout[i+1]-iat_Pout[i];
         max_nnzr = (max_nnzr>len) ?  max_nnzr:len;
      }
      printf("Final Pattern max nnzr:                      %10d\n",max_nnzr);
   }
   if (verb >= 2){
      double nt_P_d = static_cast<double>(iat_Pout[nn_A]);
      printf("Final Pattern avg nnzr:                      %10.2f\n",nt_P_d / nn_A_d);
   }

   /////////////////////////////////////////////////////
   //ofile = fopen("mat_Patt.csr","w"); if (!ofile) exit(1);
   //for (int i = 0; i < nn_A; i++){
   //    for (int j = iat_P[i]; j < iat_P[i+1]; j++){
   //       fprintf(ofile,"%10d %10d %1d\n",i+1,ja_P[j]+1,1);
   //    }
   //}
   //fclose(ofile);
   /////////////////////////////////////////////////////

   // Delete scratch
   delete [] iat_F;
   delete [] ja_F;

   return 0;
}
