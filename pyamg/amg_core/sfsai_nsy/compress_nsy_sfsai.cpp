//#include <iostream>  // to use: cout,endl
#include <algorithm> // to use: fill_n

//----------------------------------------------------------------------------------------
//
// On entry iat_FL and ja_FL hold the topology for FU as well.
//
//----------------------------------------------------------------------------------------
int compress_nsy_sfsai(const int nn, const int nt_FL, const int nt_FU,
                       int * iat_FL, int * ja_FL, int * iat_FU, int * ja_FU,
                       double * coef_FL, double * coef_FU){

   // Compress FL
   int *iat_FL_new = new int [nn+1]();
   int *ja_FL_new = new int [nt_FL]();
   double *coef_new = new double [nt_FL]();
   if (iat_FL_new == nullptr || ja_FL_new == nullptr || coef_new == nullptr){
      // Allocation error
      //std::cout << "Allocation Error in compress_nsy_sfsai" << std::endl;
      return 1;
   }
   int kk = 0;
   iat_FL_new[0] = 0;
   for (int i = 0; i < nn; i++){
      for (int j = iat_FL[i]; j < iat_FL[i+1]; j++){
         if (coef_FL[j] != 0.0){
            coef_new[kk] = coef_FL[j];
            ja_FL_new[kk] = ja_FL[j];
            kk++;
         }
      }
      iat_FL_new[i+1] = kk;
   }

   // Swap pointers for the entries of FL
   for (int i = 0; i < nt_FL; i++){
      coef_FL[i] = coef_new[i];
   }
   delete [] coef_new;

   // Compress and transpose FU
   coef_new = new double [nt_FU]();
   int *ISCR = new int[nn+1]();
   if (coef_new == nullptr || ISCR == nullptr){
      // Allocation error
      //std::cout << "Allocation Error in compress_nsy_sfsai" << std::endl;
      return 1;
   }

   // Initialize pointers
   std::fill_n(iat_FU,nn+1,0);

   // Count non-zeroes for each column of the input matrix
   for ( int i = 0; i < nn; i++ ){
      for ( int j = iat_FL[i]; j < iat_FL[i+1]; j++ ){
         if (coef_FU[j] != 0.0) iat_FU[ja_FL[j]]++;
      }
   }

   // Set pointers
   ISCR[0] = 0;
   for ( int i = 1; i < nn+1; i++ ) ISCR[i] = ISCR[i-1] + iat_FU[i-1];
   for ( int i = 0; i < nn+1; i++ ) iat_FU[i] = ISCR[i];

   // Transpose non-zero entries
   for ( int i = 0; i < nn; i++ ){
      for ( int j = iat_FL[i]; j < iat_FL[i+1]; j++ ){
         if (coef_FU[j] != 0.0){
            int ind  = ISCR[ja_FL[j]];
            ja_FU[ind] = i;
            coef_new[ind] = coef_FU[j];
            ISCR[ja_FL[j]] = ind+1;
         }
      }
   }

   // Swap pointers for the entries of FU
   for (int i = 0; i < nt_FL; i++){
      coef_FU[i] = coef_new[i];
   }
   delete [] coef_new;

   // Swap pointers for the indices of FL
   for (int i = 0; i < nn+1; i++){
      iat_FL[i] = iat_FL_new[i];
   }
   for (int i = 0; i < nt_FL; i++){
      ja_FL[i] = ja_FL_new[i];
   }
   delete [] iat_FL_new;
   delete [] ja_FL_new;

   // Delete scratch
   delete [] ISCR;

   return 0;

}
