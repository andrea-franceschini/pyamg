///////////////////////////////////////
//#include <iostream>  // to use: cout,endl
///////////////////////////////////////
#include <algorithm> // to use: fill_n

#include "iheapsort.h"

int power_patt(const int kpow, const int nnzr_max, const int nn, const int *const iat,
               const int *const ja, int &mmax, int * iat_P, const int ja_Pout_size,
               int * ja_Pout){

   int ierr = 0;

   // Allocate result arrays with a tentative size for ja_P
   int nt_P_max = std::min(6,kpow)*iat[nn];
   int * ja_P = new int [nt_P_max]();
   if (ja_P == nullptr) return 1;

   // Allocate scratch
   int *IWN = new int [nn]();
   if (IWN == nullptr) return 1;
   std::fill_n(IWN,nn,0);

   // Initialization
   mmax = 1;
   int ind = 0;
   iat_P[0] = 0;

   // Loop over rows
   int rend = iat[0];
   for (int irow = 0; irow < nn; irow++){
      // Copy the pattern of row irow of ja in ja_P
      int rstrt = rend;
      rend  = iat[irow+1];
      int ind_S = ind;
      for (int j = rstrt; j < rend; j++){
          ja_P[ind] = ja[j];
          ind++;
      }
      int ind1 = ind_S;
      // Mark corresponding entries in IWN
      for (int j = ind_S; j < ind; j++) IWN[ja_P[j]] = 1;
      // Execute all the powers
      for (int ipow = 2; ipow <= kpow; ipow++){
         // Perform product with row irow
         int ind2 = ind;
         for (int j = ind1; j < ind2; j++){
            // Compare with row jrow
            int jrow = ja_P[j];
            for (int k = iat[jrow]; k < iat[jrow+1]; k++){
               // Check if the term jcol is to be added
               int jcol = ja[k];
               if (IWN[jcol] == 0){
                  IWN[jcol] = 1;
                  ja_P[ind] = jcol;
                  ind++;
               }
            }
         }
         //------------------------- NEW CHECK ON NNZR -------------------------------
         // Check that the number of non-zeros for this row are not too many
         if (ind-ind_S > nnzr_max){
            // Sparse nullify IWN of the last added level
            for (int j = ind2; j < ind; j++) IWN[ja_P[j]] = 0;
            ind = ind2;
            break;
         }
         //---------------------------------------------------------------------------
         // Check the available memory
         if (ind + nn > nt_P_max){
            // Compute expansion factor
            double exp_fac = 1.2*static_cast<double>(nn)/static_cast<double>(irow);
            nt_P_max = static_cast<int>(exp_fac*static_cast<double>(nt_P_max));
            //////////////////////////////////////////////////////
            //std::cout << "The available memory for the pattern has been reached at row: " << irow
            //     << std::endl;
            //std::cout << "Expand the memory allocated by the factor: " << exp_fac << std::endl;
            //std::cout << "NEW number of entries allocated for the pattern: " << nt_P_max << std::endl;
            //////////////////////////////////////////////////////
            // Allocate new room
            int *ja_P_new = new int [nt_P_max]();
            if (ja_P_new == nullptr) return 1;
            // Copy old entries
            for (int j = 0; j < ind; j++) ja_P_new[j] = ja_P[j];
            // Delete old array
            delete [] ja_P;
            ja_P = ja_P_new;
         }
         ind1 = ind2;
      }
      // Sort the i-th row of A^kappa
      iheapsort(&(ja_P[ind_S]),ind-ind_S);
      // Sparse nullify IWN
      for (int j = ind_S; j < ind; j++) IWN[ja_P[j]] = 0;
      // Find the position of the diagonal term
      ind1 = ind_S;
      while (ja_P[ind1] < irow) ind1++;
      mmax = std::max(mmax,ind1-ind_S+1);
      ind = ind1 + 1;
      iat_P[irow+1] = ind;
   }

   // Delete scratch
   delete [] IWN;

   if (ind > ja_Pout_size)
   {
      ierr = -ind;
   }
   else
   {
      // Copy old entries
      for (int j = 0; j < ind; j++) ja_Pout[j] = ja_P[j];
   }
   // Delete old array
   delete [] ja_P;

   return ierr;

}
