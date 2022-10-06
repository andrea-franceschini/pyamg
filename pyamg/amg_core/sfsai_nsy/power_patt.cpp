///////////////////////////////////////
#include <iostream>  // to use: cout,endl
using std::cout;
using std::endl;
//#include <iomanip>  // to use: setprecision
//#include "wrCSRmat.h"
///////////////////////////////////////
#include <algorithm> // to use: fill_n
using std::fill_n;
using std::max;
using std::min;
#include "omp.h"

#include "iheapsort.h"

int power_patt(const int np, const int kpow, const int nnzr_max, const int nn,
               const int *const iat, const int *const ja, int &mmax,
               int *& iat_P, int *& ja_P){

   // Init error code
   int ierr = 0;

   // Allocate global scratch
   int *ridv = (int*) malloc((np+1) * sizeof(int));
   if (ridv == nullptr) return ierr = 1;

   // Init max number of non-zeroes
   mmax = 0;
   #pragma omp parallel num_threads(np) reduction(max:mmax)
   {
      int mmax_loc,ind;

      // Create the row partition
      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int mynrows,firstrow;
      if (mythid <= resto) {
         mynrows = bsize+1;
         firstrow = mythid*mynrows;
         if (mythid == resto) mynrows--;
      } else {
         mynrows = bsize;
         firstrow = mythid*bsize + resto;
      }

      // Allocate local scratches with a tentative size for ja_P
      int nt_loc = iat[firstrow+mynrows] - iat[firstrow];
      int my_nt_P_max = min(6,kpow)*nt_loc;
      int *my_iat_P = (int*) malloc((mynrows+1) * sizeof(int));
      int *my_ja_P = (int*) malloc(my_nt_P_max * sizeof(int));
      int *IWN = (int*) malloc(nn * sizeof(int));
      if (my_iat_P == nullptr || my_ja_P == nullptr || IWN == nullptr){
         #pragma omp atomic update
         ierr += 1;
      }
      #pragma omp barrier
      if (ierr != 0) goto mid_point_check;

      // Initialization
      mmax_loc = 1;
      ind = 0;
      my_iat_P[0] = 0;
      fill_n(IWN,nn,0);

      // Loop over a stripe of rows
      int rend, rstrt;
      rend = iat[firstrow];
      for (int irow = 0; irow < mynrows; irow++){
         // Copy the pattern of row irow of ja in ja_P
         rstrt = rend;
         rend  = iat[firstrow+irow+1];
         int ind_S = ind;
         for (int j = rstrt; j < rend; j++){
             my_ja_P[ind] = ja[j];
             ind++;
         }
         int ind1 = ind_S;
         // Mark corresponding entries in IWN
         for (int j = ind_S; j < ind; j++) IWN[my_ja_P[j]] = 1;
         // Execute all the powers
         for (int ipow = 2; ipow <= kpow; ipow++){
            // Perform product with row irow
            int ind2 = ind;
            for (int j = ind1; j < ind2; j++){
               // Compare with row jrow
               int jrow = my_ja_P[j];
               for (int k = iat[jrow]; k < iat[jrow+1]; k++){
                  // Check if the term jcol is to be added
                  int jcol = ja[k];
                  if (IWN[jcol] == 0){
                     IWN[jcol] = 1;
                     my_ja_P[ind] = jcol;
                     ind++;
                  }
               }
            }
            //------------------------- NEW CHECK ON NNZR -------------------------------
            // Check that the number of non-zeros for this row are not too many
            if (ind-ind_S > nnzr_max){
               // Sparse nullify IWN of the last added level
               for (int j = ind2; j < ind; j++) IWN[my_ja_P[j]] = 0;
               ind = ind2;
               break;
            }
            //---------------------------------------------------------------------------
            // Check the available memory
            if (ind + nn > my_nt_P_max){
               // Compute expansion factor
               double exp_fac = 1.2*static_cast<double>(mynrows)/static_cast<double>(irow);
               my_nt_P_max = static_cast<int>(exp_fac*static_cast<double>(my_nt_P_max));
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               //cout << "The available memory for the pattern has been reached at row: " << irow
               //     << endl;
               //cout << "Expand the memory allocated by the factor: " << exp_fac << endl;
               //cout << "NEW number of entries allocated for the pattern: " <<
               //        my_nt_P_max << endl;
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               // Allocate new room
               int *my_ja_P_new = (int*) malloc(my_nt_P_max * sizeof(int));
               if (my_ja_P_new == nullptr){
                  #pragma omp atomic update
                  ierr += 1;
                  goto mid_point_check;
               }
               // Copy old entries
               for (int j = 0; j < ind; j++) my_ja_P_new[j] = my_ja_P[j];
               // Delete old array
               free(my_ja_P);
               my_ja_P = my_ja_P_new;
            }
            ind1 = ind2;
         }
         // Sort the i-th row of A^kappa
         iheapsort(&(my_ja_P[ind_S]),ind-ind_S);
         // Sparse nullify IWN
         for (int j = ind_S; j < ind; j++) IWN[my_ja_P[j]] = 0;
         // Find the position of the diagonal term
         ind1 = ind_S;
         while (my_ja_P[ind1] < firstrow+irow) ind1++;
         mmax_loc = max(mmax_loc,ind1-ind_S+1);
         ind = ind1 + 1;
         my_iat_P[irow+1] = ind;
      }

      // Store local number of non-zeroes
      ridv[mythid+1] = ind;

      mid_point_check:

      // Delete scratch
      free(IWN);

      #pragma omp barrier
      if (ierr != 0) goto exit_point_omp;

      #pragma omp single
      {
         ridv[0] = 0;
         for (int i = 0; i < np; i++) ridv[i+1] += ridv[i];
         int nt_P = ridv[np];
         // Allocate room for filtered A
         iat_P =  (int*) malloc((nn+1) * sizeof(int));
         ja_P  =  (int*) malloc(nt_P * sizeof(int));
         if (iat_P == nullptr || ja_P == nullptr) ierr++;
      }

      // Check there is no allocation error
      if (ierr != 0) goto exit_point_omp;

      // Each thread copies its stripe of P
      ind = ridv[mythid];
      for (int i = 0; i < mynrows; i++){
         iat_P[firstrow+i] = ind;
         for (int j = my_iat_P[i]; j < my_iat_P[i+1]; j++){
            ja_P[ind] = my_ja_P[j];
            ind++;
         }
      }
      if (mythid == np-1) iat_P[nn] = ind;

      exit_point_omp:

      // Free scratches
      free(my_iat_P);
      free(my_ja_P);

      // Get max
      mmax = mmax_loc;

   }

   // Free global scratch
   free(ridv);

   return ierr;

}
