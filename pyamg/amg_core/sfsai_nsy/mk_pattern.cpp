#include <iostream>  // to use: cout,endl
#include <math.h>    // to use: fabs
#include "omp.h"
//@@@@@@@@@@@@@@@@@@@@@@@@@@
//#include <iostream>
//using std::endl;
//using std::cout;
//@@@@@@@@@@@@@@@@@@@@@@@@@@

#include "transp.h"
#include "merge_row_patt.h"
#include "power_patt.h"

int mk_pattern(const int verb, const int np, const int kpow, const double tau_pref,
               const int nnzr_max, const int nn_A, const int nt_A, const int *const iat_A,
               const int *const ja_A, const double *const coef_A, int &mmax,
               int *&iat_P, int *&ja_P){

   // Init error code
   int ierr = 0;

   // Extract diagonals from A
   double *diag_A = (double*) malloc(nn_A*sizeof(double));
   if (diag_A == nullptr) return 1;
   #pragma omp parallel for num_threads(np)
   for (int i = 0; i < nn_A; i++){
      int ind = iat_A[i];
      int jcol = ja_A[ind];
      //@@@@ HERE WE ASSUME THAT EACH ROW CONTAINS A NON-ZERO DIAGONAL ENTRY
      while (jcol < i){
         ind++;
         jcol = ja_A[ind];
      }
      diag_A[i] = fabs(coef_A[ind]);
   }

   // Allocate scratch
   int *ridv =  (int*) malloc((np+1) * sizeof(int));
   if (ridv == nullptr) return 1;

   // Filter A
   int nt_FA;
   int *iat_FA;
   int *ja_FA;
   #pragma omp parallel num_threads(np)
   {
      // Create the row partition
      int mythid = omp_get_thread_num();
      int bsize = nn_A/np;
      int resto = nn_A%np;
      int mynrows,firstrow;
      if (mythid <= resto) {
         mynrows = bsize+1;
         firstrow = mythid*mynrows;
         if (mythid == resto) mynrows--;
      } else {
         mynrows = bsize;
         firstrow = mythid*bsize + resto;
      }

      // Allocate local scratches
      int nnz_prev = iat_A[firstrow+mynrows] - iat_A[firstrow];
      int *my_iat_FA = (int*) malloc((mynrows+1) * sizeof(int));
      int *my_ja_FA  = (int*) malloc(nnz_prev * sizeof(int));
      if (my_iat_FA == nullptr || my_ja_FA == nullptr){
         #pragma omp atomic update
         ierr += 1;
      }
      #pragma omp barrier
      if (ierr != 0) goto exit_point_omp_1;

      // Loop on each stripe of rows
      int nt_FA_loc, k;
      nt_FA_loc = 0;
      k = 0;
      for (int i = firstrow; i < firstrow+mynrows; i++){
         my_iat_FA[k] = nt_FA_loc;
         k++;
         double fac = tau_pref*diag_A[i];
         for (int j = iat_A[i]; j < iat_A[i+1]; j++){
            int jcol = ja_A[j];
            if (pow(fabs(coef_A[j]),2) > fac*diag_A[jcol]){
               my_ja_FA[nt_FA_loc] = jcol;
               nt_FA_loc++;
            }
         }
      }
      my_iat_FA[k] = nt_FA_loc;
      ridv[mythid+1] = nt_FA_loc;

      // Synchronize threads
      #pragma omp barrier

      #pragma omp single
      {
         ridv[0] = 0;
         for (int i = 0; i < np; i++) ridv[i+1] += ridv[i];
         nt_FA = ridv[np];
         // Allocate room for filtered A
         iat_FA =  (int*) malloc((nn_A+1) * sizeof(int));
         ja_FA  =  (int*) malloc(nt_FA * sizeof(int));
         if (iat_FA == nullptr || ja_FA == nullptr) ierr = 1;
      }

      // Check there is no allocation error
      if (ierr != 0) goto exit_point_omp_1;

      // Each thread copies its stripe of filtered A
      int ind;
      ind = ridv[mythid];
      for (int i = 0; i < mynrows; i++){
         iat_FA[firstrow+i] = ind;
         for (int j = my_iat_FA[i]; j < my_iat_FA[i+1]; j++){
            ja_FA[ind] = my_ja_FA[j];
            ind++;
         }
      }
      if (mythid == np-1) iat_FA[nn_A] = ind;

      exit_point_omp_1:

      // Free local scratch
      free(my_iat_FA);
      free(my_ja_FA);

   }

   // Check for error
   if (ierr != 0) return ierr = 1;

   // Free scratches
   free(diag_A);

   double nn_A_d = static_cast<double>(nn_A);
   double nt_FA_d = static_cast<double>(nt_FA);
   if (verb >= 1){
      fprintf(stdout,"Initial Filtered Pattern avg nnzr:           %10.2f\n",
              nt_FA_d / nn_A_d);
   }

   /////////////////////////////////////////////////////
   //FILE *ofile = fopen("mat_InitPatt.csr","w"); if (!ofile) exit(1);
   //for (int i = 0; i < nn_A; i++){
   //    for (int j = iat_FA[i]; j < iat_FA[i+1]; j++){
   //       fprintf(ofile,"%10d %10d %1d\n",i+1,ja_FA[j]+1,1);
   //    }
   //}
   //fflush(ofile);
   //fclose(ofile);
   /////////////////////////////////////////////////////

   // Transpose filtered pattern
   int *iat_FT;
   int *ja_FT;
   ierr = transp_patt(np,nn_A,nn_A,iat_FA,ja_FA,iat_FT,ja_FT);
   if (ierr != 0) return ierr = 1;

   // Merge filtered and filtered transposed patterns
   int nt_F;
   int *iat_F;
   int *ja_F;
   #pragma omp parallel num_threads(np)
   {
      // Create the row partition
      int mythid = omp_get_thread_num();
      int bsize = nn_A/np;
      int resto = nn_A%np;
      int mynrows,firstrow;
      if (mythid <= resto) {
         mynrows = bsize+1;
         firstrow = mythid*mynrows;
         if (mythid == resto) mynrows--;
      } else {
         mynrows = bsize;
         firstrow = mythid*bsize + resto;
      }

      // Allocate local scratches
      int nt_F_loc;
      int max_nt_F_loc = iat_FA[firstrow+mynrows] - iat_FA[firstrow] +
                         iat_FT[firstrow+mynrows] - iat_FT[firstrow];
      int *my_iat_F = (int*) malloc((mynrows+1) * sizeof(int));
      int *my_ja_F = (int*) malloc(max_nt_F_loc * sizeof(int));
      if (my_iat_F == nullptr || my_ja_F == nullptr){
         #pragma omp atomic update
         ierr += 1;
      }
      #pragma omp barrier
      if (ierr != 0) goto exit_point_omp_2;

      int ind, iend_FA, iend_FT;
      ind = 0;
      iend_FA = iat_FA[firstrow];
      iend_FT = iat_FT[firstrow];
      for (int i = 0; i < mynrows; i++){
         int len_out;

         // Set pointer to output
         my_iat_F[i] = ind;
         // Set pointer to input 1
         int istrt_FA = iend_FA;
         iend_FA = iat_FA[firstrow+i+1];
         int len_FA = iend_FA-istrt_FA;
         // Set pointer to input 2
         int istrt_FT = iend_FT;
         iend_FT = iat_FT[firstrow+i+1];
         int len_FT = iend_FT-istrt_FT;
         // Merge the rows
         merge_row_patt(len_FA,&(ja_FA[istrt_FA]),len_FT,&(ja_FT[istrt_FT]),
                        len_out,&(my_ja_F[ind]));
         // Update pointer
         ind += len_out;
      }
      my_iat_F[mynrows] = ind;
      nt_F_loc = ind;
      ridv[mythid+1] = nt_F_loc;

      // Synchronize threads
      #pragma omp barrier

      #pragma omp single
      {
         ridv[0] = 0;
         for (int i = 0; i < np; i++) ridv[i+1] += ridv[i];
         nt_F = ridv[np];
         // Allocate room for F
         iat_F =  (int*) malloc((nn_A+1) * sizeof(int));
         ja_F  =  (int*) malloc(nt_F * sizeof(int));
         if (iat_F == nullptr || ja_F == nullptr) ierr = 1;
      }

      // Check there is no allocation error
      if (ierr != 0) goto exit_point_omp_2;

      // Each thread copies its stripe of F
      ind = ridv[mythid];
      for (int i = 0; i < mynrows; i++){
         iat_F[firstrow+i] = ind;
         for (int j = my_iat_F[i]; j < my_iat_F[i+1]; j++){
            ja_F[ind] = my_ja_F[j];
            ind++;
         }
      }
      if (mythid == np-1) iat_F[nn_A] = ind;

      exit_point_omp_2:

      // Free scratches
      free(my_iat_F);
      free(my_ja_F);

   }
   if (ierr != 0) return ierr = 1;

   // Delete intermediate scratches
   free(ridv);
   free(iat_FT);
   free(ja_FT);
   free(iat_FA);
   free(ja_FA);

   double nt_F_d = static_cast<double>(nt_F);
   if (verb >= 1){
      fprintf(stdout,"Final Filtered Symmetrized Pattern avg nnzr: %10.2f\n",
              nt_F_d / nn_A_d);
   }

   /////////////////////////////////////////////////////
   //{
   //FILE *ofile = fopen("mat_FiltPatt.csr","w"); if (!ofile) exit(1);
   //for (int i = 0; i < nn_A; i++){
   //    for (int j = iat_F[i]; j < iat_F[i+1]; j++){
   //       fprintf(ofile,"%10d %10d %1d\n",i+1,ja_F[j]+1,1);
   //    }
   //}
   //fflush(ofile);
   //fclose(ofile);
   //}
   /////////////////////////////////////////////////////

   // Compute the power of the final filtered pattern
   ierr = power_patt(np,kpow,nnzr_max,nn_A,iat_F,ja_F,mmax,iat_P,ja_P);
   if (ierr != 0) return 1;

   double nt_P_d = static_cast<double>(iat_P[nn_A]);
   if (verb >= 1){
      fprintf(stdout,"Final Pattern avg nnzr:                      %10.2f\n",
              nt_P_d / nn_A_d);
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
   free(iat_F);
   free(ja_F);

   return 0;
}
