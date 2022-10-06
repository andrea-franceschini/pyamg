#include <algorithm> // to use: fill_n
#include "omp.h"

//----------------------------------------------------------------------------------------
//
// On entry iat_FL and ja_FL hold the topology for FU as well.
//
//----------------------------------------------------------------------------------------
int compress_nsy_sfsai(const int np, const int nn, const int *pt_FL, const int *pt_FU,
                       int *&iat_FL, int *&ja_FL, int *&iat_FU, int *&ja_FU,
                       double *&coef_FL, double *&coef_FU){

   // Compress FL and FU
   int *iat_FL_new = (int*) malloc((nn+1) * sizeof(int));
   int *ja_FL_new = (int*) malloc(pt_FL[np] * sizeof(int));
   double *coef_FL_new = (double*) malloc(pt_FL[np] * sizeof(double));
   iat_FU = (int*) malloc((nn+1) * sizeof(int));
   ja_FU = (int*) malloc(pt_FU[np] * sizeof(int));
   double *coef_FU_new = (double*) malloc(pt_FU[np] * sizeof(double));
   if (iat_FL_new == nullptr || ja_FL_new == nullptr || coef_FL_new == nullptr ||
       iat_FU == nullptr || ja_FU == nullptr || coef_FU_new == nullptr) return 1;

   #pragma omp parallel num_threads(np)
   {
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

      int kk = pt_FL[mythid];
      int ll = pt_FU[mythid];
      for (int i = 0; i < mynrows; i++){
         iat_FL_new[firstrow+i] = kk;
         iat_FU[firstrow+i] = ll;
         for (int j = iat_FL[firstrow+i]; j < iat_FL[firstrow+i+1]; j++){
            if (coef_FL[j] != 0.0){
               coef_FL_new[kk] = coef_FL[j];
               ja_FL_new[kk] = ja_FL[j];
               kk++;
            }
            if (coef_FU[j] != 0.0){
               coef_FU_new[ll] = coef_FU[j];
               ja_FU[ll] = ja_FL[j];
               ll++;
            }
         }
      }
      if (mythid == np-1){
         iat_FL_new[nn] = kk;
         iat_FU[nn] = ll;
      }

   }

   // Swap pointers for the entries of FL
   free(coef_FL);
   coef_FL = coef_FL_new;

   // Swap pointers for the entries of FU
   free(coef_FU);
   coef_FU = coef_FU_new;

   // Swap pointers for topology of FL
   free(iat_FL);
   free(ja_FL);
   iat_FL = iat_FL_new;
   ja_FL = ja_FL_new;

   return 0;

}
