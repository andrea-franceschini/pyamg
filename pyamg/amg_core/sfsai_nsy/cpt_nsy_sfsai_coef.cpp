#define DEBUG false
//#define DEBUG true
#include <cstring>
///////////////////////////////////////////

#include <iostream>
using std::cout;
using std::endl;

#include "omp.h"
#include <lapacke.h>  // to use: dgetrf,dgetrs

#include "inl_blas1.h"
#include "gather_fullsys.h"
#include "fsai_ScalFact.h"

int cpt_nsy_sfsai_coef(const int np, const double tau, const int mmax, const int nn,
                       const int *const iat_A, const int *const ja_A,
                       const double *const coef_A, const int *const iat_F,
                       const int *const ja_F, int *nt_FL, int *nt_FU,
                       double *coef_FL, double *coef_FU){

   // Init error code
   int ierr = 0;

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

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      FILE *ofile;
      if (DEBUG){
         char nome_out[30];
         char label[4];
         strcpy(nome_out,"Log_cpt_coef_");
         sprintf(label, "%2d",mythid);
         strcat(nome_out,label);
         ofile = fopen(nome_out,"w"); if (!ofile) exit(1);
      }
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      // Allocate local scratches
      double *full_A = (double*) malloc((mmax-1)*(mmax-1) * sizeof(double));
      double *rhs_L = (double*) malloc(mmax * sizeof(double));
      double *rhs_U = (double*) malloc(mmax * sizeof(double));
      lapack_int *ipvt = (lapack_int*) malloc((mmax-1) * sizeof(lapack_int));
      double *full_A_sav = (double*) malloc((mmax-1)*(mmax-1) * sizeof(double));
      double *rhs_L_sav = (double*) malloc(mmax * sizeof(double));
      double *rhs_U_sav = (double*) malloc(mmax * sizeof(double));
      if (full_A == nullptr || rhs_L == nullptr || rhs_U == nullptr ||
          ipvt == nullptr || full_A_sav == nullptr || rhs_L_sav == nullptr ||
          rhs_U_sav == nullptr){
          #pragma omp atomic update
          ierr++;
          goto exit_point_omp;
      }

      // Init the number of non-zeros of this stripe of FL and FU
      int my_nt_FL, my_nt_FU;
      my_nt_FL = iat_F[firstrow+mynrows] - iat_F[firstrow];
      my_nt_FU = iat_F[firstrow+mynrows] - iat_F[firstrow];

      // Loop over the rows of F
      int iend_F;
      iend_F = iat_F[firstrow];
      for (int irow = 0; irow < mynrows; irow++){

          // Update pointers
          int istart_F = iend_F;
          iend_F = iat_F[firstrow+irow+1];

          // Get current number of rows
          int mrow = (iend_F - istart_F) - 1;

          // Compute coefficients
          double diag_entry;
          if (mrow == 0){

             // Just retrieve diagonal entry
             int jj = iat_A[firstrow+irow];
             diag_entry = 0.0;
             while (ja_A[jj] < firstrow+irow){
                jj++;
                if (jj == iat_A[firstrow+irow+1]){
                   goto while_exit;
                }
             }
             if (ja_A[jj] == firstrow+irow) diag_entry = coef_A[jj];
             while_exit:;

          } else {

             // Gather the coefficients of the full local systems
             bool null_L;
             bool null_U;
             gather_fullsys(mrow,&(ja_F[istart_F]),nn,iat_A,ja_A,coef_A,diag_entry,full_A,
                            rhs_L,rhs_U,null_L,null_U);
             //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
             if (DEBUG){
                fprintf(ofile,"-------------------------\n");
                fprintf(ofile,"irow = %d\n\n",irow);
                fprintf(ofile,"mrow = %d\n\n",mrow);
                fprintf(ofile,"vecinc = ");
                for (int i = 0; i < mrow; i++) fprintf(ofile,"%d ",ja_F[istart_F+i]);
                fprintf(ofile,"\n\n");
                fprintf(ofile,"Diag_A = %10e\n\n",diag_entry);
                fprintf(ofile,"A:\n");
                int k = 0;
                for (int i = 0; i < mrow; i++){
                   for (int j = 0; j < mrow; j++) fprintf(ofile,"%10e ",full_A[k++]);
                   fprintf(ofile,"\n");
                }
                fprintf(ofile,"\n");
                fprintf(ofile,"rhs_U:\n");
                for (int j = 0; j < mrow; j++) fprintf(ofile,"%10e ",rhs_U[j]);
                fprintf(ofile,"\n\n");
                fprintf(ofile,"rhs_L:\n");
                for (int j = 0; j < mrow; j++) fprintf(ofile,"%10e ",rhs_L[j]);
                fprintf(ofile,"\n\n");
             }
             //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

             // Backup system and rhs
             for (int k = 0; k < mrow*mrow; k++) full_A_sav[k] = full_A[k];
             for (int k = 0; k < mrow; k++) rhs_L_sav[k] = rhs_L[k];
             for (int k = 0; k < mrow; k++) rhs_U_sav[k] = rhs_U[k];

             // Factorize the dense matrix
             if (!null_L || !null_U){
                lapack_int ierr_L = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,mrow,mrow,full_A,mrow,
                                                   ipvt);
                if (ierr_L != 0){
                   // Lapack error
                   printf("DGETRF: error %d at row %d for thread %d\n",
                          ierr_L,firstrow+irow,mythid);
                   #pragma omp atomic update
                   ierr++;
                   goto exit_point_omp;
                }
             }
             
             // Compute coefficients of the irow-th row of FL / column of FU
             if (!null_L){
                lapack_int ierr_L = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'T',mrow,1,full_A,mrow,
                                                   ipvt,rhs_L,1);
                if (ierr_L != 0){
                   // Lapack error
                   printf("DGETRS: error %d at row %d for thread %d\n",
                          ierr_L,firstrow+irow,mythid);
                   #pragma omp atomic update
                   ierr++;
                   goto exit_point_omp;
                }
             }
             if (!null_U){
                lapack_int ierr_L = LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',mrow,1,full_A,mrow,
                                                 ipvt,rhs_U,1);
                if (ierr_L != 0){
                   // Lapack error
                   printf("DGETRS: error %d at row %d for thread %d\n",
                          ierr_L,firstrow+irow,mythid);
                   #pragma omp atomic update
                   ierr++;
                   goto exit_point_omp;
                }
             }
             //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
             if (DEBUG){
                fprintf(ofile,"sol_U:\n");
                for (int j = 0; j < mrow; j++) fprintf(ofile,"%10e ",rhs_U[j]);
                fprintf(ofile,"\n\n");
                fprintf(ofile,"sol_L:\n");
                for (int j = 0; j < mrow; j++) fprintf(ofile,"%10e ",rhs_L[j]);
                fprintf(ofile,"\n\n");
             }
             //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          }

          // Filter small entries
          rhs_L[mrow] = 1.0;
          double absTol = tau*inl_dnrm2(mrow+1,rhs_L,1);
          for (int k = 0; k < mrow; k++){
             if (fabs(rhs_L[k]) <= absTol){
                rhs_L[k] = 0.0;
                my_nt_FL--;
             }
          }
          rhs_U[mrow] = 1.0;
          absTol = tau*inl_dnrm2(mrow+1,rhs_U,1);
          for (int k = 0; k < mrow; k++){
             if (fabs(rhs_U[k]) <= absTol){
                rhs_U[k] = 0.0;
                my_nt_FU--;
             }
          }

          // Compute scaling factor
          double scal_fac = fsai_ScalFact(mrow,diag_entry,rhs_L_sav,rhs_U_sav,full_A_sav,
                                          rhs_L,rhs_U,full_A);
          //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          if (DEBUG) fprintf(ofile,"SCAL_FAC: %e\n\n",scal_fac);
          //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          // Check zero diagonal
          double check_val = fabs(scal_fac / diag_entry);
          if ( check_val < 1.0e-10 ){
             cout << "SMALL DIAGONAL = " << check_val << " IN ROW: " << irow << endl;
          } 

          // Scale and store lower part
          double fac = 1.0 / sqrt(fabs(scal_fac));
          //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          if (DEBUG) fprintf(ofile,"INV_SCAL_FAC: %e\n\n",fac);
          //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          for (int k = 0; k < mrow+1; k++) coef_FL[istart_F+k] = rhs_L[k]*fac;

          // Scale and store upper part
          if (scal_fac < 0.0) fac = -fac;
          for (int k = 0; k < mrow+1; k++) coef_FU[istart_F+k] = rhs_U[k]*fac;

      } // End row loop

      // Save current number of non-zeroes
      nt_FL[mythid] = my_nt_FL;
      nt_FU[mythid] = my_nt_FU;

      exit_point_omp:

      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if (DEBUG) fclose(ofile);
      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      // Deallocate local scratches
      free(rhs_U_sav);
      free(rhs_L_sav);
      free(full_A_sav);
      free(ipvt);
      free(rhs_U);
      free(rhs_L);
      free(full_A);

   }

   return ierr;

}
