#include <iostream>
#include <algorithm> // to use: fill_n

int transp_Patt(const int nrows, const int ncols, const int *const iat, const int *const ja,
                int *&iat_T, int*&ja_T){

   // Allocate output and scracth
   iat_T = new int [ncols+1]();
   int nterm = iat[nrows];
   ja_T = new int [nterm]();
   int *ISCR = new int[ncols+1]();
   if (iat_T == nullptr || ja_T == nullptr || ISCR == nullptr){
      // Allocation error
      //std::cout << "Allocation Error in transp_Patt" << std::endl;
      return 1;
   }

   // Initialize pointers
   std::fill_n(iat_T,ncols+1,0);

   // Count non-zeroes for each column of the input matrix
   for ( int i = 0; i < nrows; i++ ){
      for ( int j = iat[i]; j < iat[i+1]; j++ ) iat_T[ja[j]]++;
   }

   // Set pointers
   ISCR[0] = 0;
   for ( int i = 1; i < ncols+1; i++ ) ISCR[i] = ISCR[i-1] + iat_T[i-1];
   for ( int i = 0; i < ncols+1; i++ ) iat_T[i] = ISCR[i];

   // Transpose column indices
   for ( int i = 0; i < nrows; i++ ){
      for ( int j = iat[i]; j < iat[i+1]; j++ ){
         int ind  = ISCR[ja[j]];
         ja_T[ind] = i;
         ISCR[ja[j]] = ind+1;
      }
   }

   // Deallocate scratch
   delete [] ISCR;

   return 0;

}
