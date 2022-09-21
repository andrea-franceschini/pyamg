#include <stdlib.h>
#include <omp.h>
#include <limits>
#include "cblas.h"    // to use: DGEMV
#include "lapacke.h"  // to use: DGEQRF and DORGQR
//////////////////////////////////
#include "inl_blas1.h"
#include <iostream>
#include <stdio.h>
//////////////////////////////////

#include "DebEnv.h"
#include "EMIN_parm.h"

#define PRINT_LOC_INFO 0
#define PRINT_LOC_INFO_ALL 0

/*****************************************************************************************
 *
 * This function gathers from the prolongation pattern (given by rows) and the test vector
 * array TV, the block diagonal matrix B that is used to enforce the constraint. B is
 * immediately factorized with QR and only Q is returned.
 *
 * Error code:
 *
 * 0 ---> successful run
 * 1 ---> allocation error for global scratches
 * 2 ---> allocation error for private scratches
 * 3 ---> allocation error for the final output
 * 4 ---> lapack error
 *
*****************************************************************************************/
int gather_B_QR(const int np, const double condmax, const int nn, const int nn_C,
                const int ntv, const int *fcnode, const int *iat_patt, const int *ja_patt,
                const double *const *TV, double *&mat_Q, double *coef_P0)
{
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// g E tau POSSONO ESSERE DEI VETTORI LOCALI DI DIMENSIONE NTV CHE POI VENGONO CANCELLATI
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //FILE *bbf = fopen("LOG_gather","w");
   //FILE *of = fopen("P0_prima","w");
   //for (int i = 0; i < iat_patt[nn]; i++) fprintf(of,"%20.11e\n",coef_P0[i]);
   //fclose(of);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Init error code
   int ierr = 0;

   // Allocate Q
   int nrows_Q = iat_patt[nn];
   int nterm_Q = ntv*nrows_Q;
   mat_Q = (double*) malloc( nterm_Q*sizeof(double) );
   if (mat_Q == nullptr) return ierr = 3;

   // Allocate shared scratches
   int *c2glo = (int*) malloc( nn_C*sizeof(int) );
   double *vec_g = (double*) malloc( nrows_Q*sizeof(double) );
   if (c2glo == nullptr || vec_g == nullptr)  return ierr = 1;

   #pragma omp parallel num_threads(np)
   {
      // Get thread ID and column partition
      int mythid = omp_get_thread_num();
      int bsize = nn/np;
      int resto = nn%np;
      int firstcol, ncolth, lastcol;
      if (mythid <= resto) {
         ncolth = bsize+1;
         firstcol = mythid*ncolth;
         if (mythid == resto) ncolth--;
      } else {
         ncolth = bsize;
         firstcol = mythid*bsize + resto;
      }
      lastcol = firstcol + ncolth;

      // Estimate number of rows and first position for this chunk of columns
      int pos_g = iat_patt[firstcol];
      int pos_Q = pos_g*ntv;

      // Set proper position in Q and g
      double *g_scr = &(vec_g[pos_g]);
      double *BB_scr = &(mat_Q[pos_Q]);

      // Find the largest number of rows in a block
      int nrmax_blk = 0;
      int istart, iend;
      iend = iat_patt[firstcol];
      for (int i = firstcol+1; i <= lastcol; i++){
         istart = iend;
         iend = iat_patt[i];
         nrmax_blk = std::max(nrmax_blk,iend-istart);
      }

      // Query work space for DGEQRF and DORGQR
      lapack_int ierr_lapack;
      lapack_int l_nn = static_cast<lapack_int>(std::max(ntv,nrmax_blk));
      lapack_int l_mm = static_cast<lapack_int>(ntv);
      lapack_int l_kk = static_cast<lapack_int>(ntv);
      double query_work_1;
      double query_work_2;
      double query_work_3;
      double *SIGMA = nullptr;      // (DI DIMENSIONE l_nn == nr_BB_loc)
      double *dummy_U = nullptr;    // (FASULLO NON VIENE USATO)
      double *VT = nullptr;         // (DI DIMENSIONE l_mm*l_mm = ntv*ntv)
      double *tau = nullptr;
      double *work = nullptr;

      lapack_int lwork = -1;
      ierr_lapack = LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR,l_nn,l_mm,BB_scr,l_nn,tau,
                                        &query_work_1,lwork);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      ierr_lapack = LAPACKE_dorgqr_work(LAPACK_COL_MAJOR,l_nn,l_mm,l_kk,BB_scr,l_nn,tau,
                                        &query_work_2,lwork);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      ierr_lapack = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,'o','s',l_nn,l_mm,BB_scr,l_nn,
                                        SIGMA,dummy_U,l_nn,VT,l_mm,&query_work_3,lwork);
      if (ierr_lapack != 0){
         #pragma omp atomic write
         ierr = 4;
      }
      if (ierr > 0) goto exit_pragma;
      query_work_1 = std::max(query_work_1,query_work_2);
      query_work_3 = std::max(query_work_3,static_cast<double>(ntv));
      lwork = static_cast<lapack_int>(std::max(query_work_1,query_work_3));

      // Allocate private workspace
      SIGMA = (double*) malloc( std::max(ntv,nrmax_blk)*sizeof(double) );
      VT = (double*) malloc( ntv*ntv*sizeof(double) );
      tau = (double*) malloc( (ncolth*ntv)*sizeof(double) );
      work = (double*) malloc( lwork*sizeof(double) );
      if (SIGMA == nullptr || VT == nullptr || tau == nullptr || work == nullptr){
         #pragma omp atomic write
         ierr = 2;
      }
      if (ierr > 0) goto exit_pragma;

      // Create mapping from coarse node numbering to global (original) numbering
      #pragma omp for
      for (int i = 0; i < nn; i++){
         int k = fcnode[i];
         if (k >= 0) c2glo[k] = i;
      }

      // Loop over the current chunk of columns
      int ind_g, ind_BB, ind_tau, nnz_BB;
      ind_g = 0;
      ind_BB = 0;
      ind_tau = 0;
      int istart_patt, iend_patt;
      iend_patt = iat_patt[firstcol];
      for (int icol = firstcol; icol < lastcol; icol++){
         // Check that this is a FINE node
         if (fcnode[icol] < 0){

            istart_patt = iend_patt;
            iend_patt = iat_patt[icol+1];
            int nr_BB_loc = iend_patt-istart_patt;

            // Check that the row is not empty
            if (nr_BB_loc > 0){
               // Copy TV into BB_scr
               int k = ind_BB;
               for (int i = istart_patt; i < iend_patt; i++){
                  int i_F = c2glo[ja_patt[i]];
                  int kk = k;
                  for (int j = 0; j < ntv; j++){
                     BB_scr[kk] = TV[i_F][j];
                     kk += nr_BB_loc;
                  }
                  k++;
               }
               // Copy TV into g
               for (int j = 0; j < ntv; j++) g_scr[ind_g+j] = TV[icol][j];

               // Compute g = g - BB^T*coef_P0
               cblas_dgemv(CblasColMajor,CblasTrans,nr_BB_loc,ntv,-1.0,&(BB_scr[ind_BB]),
                           nr_BB_loc,&(coef_P0[istart_patt]),1,1.0,&(g_scr[ind_g]),1);
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               if (DEBUG){
                  fprintf(DebEnv.t_logfile[mythid],"-----------------------------------\n");
                  fprintf(DebEnv.t_logfile[mythid],"GLOBAL ROW OF A (ICOL): %d\n",icol);
                  fprintf(DebEnv.t_logfile[mythid],"Size of B_loc: %d %d\n",nr_BB_loc,ntv);
                  fprintf(DebEnv.t_logfile[mythid],"\nB_loc:\n");
                  for (int i = 0; i < nr_BB_loc; i++){
                     for (int j = 0; j < ntv; j++)
                        fprintf(DebEnv.t_logfile[mythid]," %17.10e",BB_scr[ind_BB+j*nr_BB_loc+i]);
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                  }
                  fprintf(DebEnv.t_logfile[mythid],"\nP0 init:\n");
                  for (int j = 0; j < nr_BB_loc; j++)
                     fprintf(DebEnv.t_logfile[mythid]," %15.6e",coef_P0[istart_patt+j]);
                  fprintf(DebEnv.t_logfile[mythid],"\nrhs:\n");
                  for (int j = 0; j < ntv; j++)
                     fprintf(DebEnv.t_logfile[mythid]," %15.6e",g_scr[ind_g+j]);
                  fprintf(DebEnv.t_logfile[mythid],"\n");
                  fflush(DebEnv.t_logfile[mythid]);
               }
               //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

               //----------------------------------+
               // Try the solution of full-rank BB |
               //----------------------------------+

               bool FAIL_QR = false;
               int i_lpk_1, i_lpk_2, i_lpk_3;
               if (nr_BB_loc >= ntv){

                  // Perform QR on BB
                  l_nn = static_cast<lapack_int>(nr_BB_loc);
                  lapack_int l_ll = l_nn;
                  i_lpk_1 = LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR,l_nn,l_mm,
                                    &(BB_scr[ind_BB]),l_ll,&(tau[ind_tau]),work,lwork);

                  // Check conditioning of the resulting R
                  double max_DR = 0.0;
                  double min_DR = std::numeric_limits<double>::max();
                  for (int kk = 0; kk < l_mm; kk++){
                     max_DR = std::max(max_DR,std::abs(BB_scr[ind_BB+kk*l_ll+kk]));
                     min_DR = std::min(min_DR,std::abs(BB_scr[ind_BB+kk*l_ll+kk]));
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
                  if (DEBUG){
                     fprintf(DebEnv.t_logfile[mythid],"DIAG R: ");
                     for (int kk = 0; kk < l_mm; kk++)
                        fprintf(DebEnv.t_logfile[mythid]," %15.6e",std::abs(BB_scr[ind_BB+kk*l_ll+kk]));
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                     fprintf(DebEnv.t_logfile[mythid],
                             "%6d max_DR %15.6e min_DR %15.6e COND_1 %15.6e\n",icol,
                             max_DR,min_DR,max_DR/min_DR);
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2
                  if (max_DR / min_DR > condmax){

                     //------------------------------------+
                     // BB is rank deficient use SVD below |
                     //------------------------------------+

                     FAIL_QR = true;
                     // Reload TV in BB_scr
                     int k = ind_BB;
                     for (int i = istart_patt; i < iend_patt; i++){
                        int i_F = c2glo[ja_patt[i]];
                        int kk = k;
                        for (int j = 0; j < ntv; j++){
                           BB_scr[kk] = TV[i_F][j];
                           kk += nr_BB_loc;
                        }
                        k++;
                     }
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     if ( PRINT_LOC_INFO ) std::cout << icol <<
                        " conditioning larger than threshold: "
                        << max_DR / min_DR << " > " << condmax << std::endl;
                     if (DEBUG) fprintf(DebEnv.t_logfile[mythid],
                            "%6d COND_LRG: max %15.6e min %15.6e cond %15.6e\n",
                            icol,max_DR,min_DR,max_DR / min_DR);
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  } else {

                     //-----------------+
                     // BB is full rank |
                     //-----------------+

                     // Solve transposed triangular system g = inv(RR')*g
                     i_lpk_2 = LAPACKE_dtrtrs_work(LAPACK_COL_MAJOR,'U','T','N',l_mm,1,
                                       &(BB_scr[ind_BB]),l_ll,&(g_scr[ind_g]),l_mm);

                     // Transform QQ from Householder rotation to standard form
                     i_lpk_3 = LAPACKE_dorgqr_work(LAPACK_COL_MAJOR,l_nn,l_mm,l_kk,
                                       &(BB_scr[ind_BB]),l_ll,&(tau[ind_tau]),work,lwork);
                     if (i_lpk_1 || i_lpk_2 || i_lpk_3){
                        #pragma omp atomic write
                        ierr = 4;
                        goto exit_loop_icol;
                     }

                     // Update coef_P0 to ensure the TV constraint: coef_P0 += QQ*g
                     cblas_dgemv(CblasColMajor,CblasNoTrans,nr_BB_loc,ntv,1.0,&(BB_scr[ind_BB]),
                                 nr_BB_loc,&(g_scr[ind_g]),1,1.0,&(coef_P0[istart_patt]),1);

                     // Compute the number of entries of Q
                     nnz_BB = nr_BB_loc*ntv;
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                     if (DEBUG){
                        fprintf(DebEnv.t_logfile[mythid],"\nQ:\n");
                        for (int i = 0; i < nr_BB_loc; i++){
                           for (int j = 0; j < std::min(ntv,nr_BB_loc); j++)
                              fprintf(DebEnv.t_logfile[mythid]," %17.10e",BB_scr[ind_BB+j*nr_BB_loc+i]);
                           fprintf(DebEnv.t_logfile[mythid],"\n");
                        }
                        fprintf(DebEnv.t_logfile[mythid],"\nDP:\n");
                        for (int j = 0; j < nr_BB_loc; j++)
                          fprintf(DebEnv.t_logfile[mythid]," %15.6e",coef_P0[istart_patt+j]);
                        fprintf(DebEnv.t_logfile[mythid],"\n");
                        fflush(DebEnv.t_logfile[mythid]);
                     }
                     //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  }

               }
               //--------------------------------------------+
               // Perform least square solution if necessary |
               //--------------------------------------------+

               if ( (nr_BB_loc < ntv) || FAIL_QR){

                  if ( PRINT_LOC_INFO )
                     std::cout << icol << " Solution for RANK deficient system " << std::endl;

                  // Compute SVD of BB
                  l_nn = static_cast<lapack_int>(nr_BB_loc);
                  lapack_int l_ll = l_nn;
                  i_lpk_1 = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR,'o','s',l_nn,l_mm,
                                                &(BB_scr[ind_BB]),l_ll,SIGMA,dummy_U,
                                                l_ll,VT,l_mm,work,lwork);
                  if (i_lpk_1){
                     #pragma omp atomic write
                     ierr = 4;
                     goto exit_loop_icol;
                  }

                  // Compute the rank of BB using SIGMA
                  int rank_BB = 1;
                  while (SIGMA[0]/SIGMA[rank_BB] < condmax){
                     rank_BB++;
                     if (rank_BB == ntv){
                        rank_BB++;
                        break;
                     }
                  }
                  rank_BB--;

                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  if (DEBUG){
                     if (nr_BB_loc < ntv) fprintf(DebEnv.t_logfile[mythid],"FAT B\n");
                     //fprintf(DebEnv.t_logfile[mythid],"\nRANK B: %d\n",rank_BB);
                     fprintf(DebEnv.t_logfile[mythid],"\n%6d COND_RID %15.6e RANK B %3d\n",icol,
                             SIGMA[0]/SIGMA[rank_BB],rank_BB);
                     fprintf(DebEnv.t_logfile[mythid],"\nU = Q:\n");
                     for (int i = 0; i < nr_BB_loc; i++){
                        for (int j = 0; j < std::min(ntv,nr_BB_loc); j++)
                           fprintf(DebEnv.t_logfile[mythid]," %15.6e",BB_scr[ind_BB+j*nr_BB_loc+i]);
                        fprintf(DebEnv.t_logfile[mythid],"\n");
                     }
                     fprintf(DebEnv.t_logfile[mythid],"\nSIGMA:\n");
                     for (int j = 0; j < std::min(ntv,nr_BB_loc); j++)
                        fprintf(DebEnv.t_logfile[mythid]," %15.6e",SIGMA[j]);
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                     fprintf(DebEnv.t_logfile[mythid],"\nVT:\n");
                     for (int i = 0; i < std::min(nr_BB_loc,ntv); i++){
                        for (int j = 0; j < ntv; j++)
                           fprintf(DebEnv.t_logfile[mythid]," %15.6e",VT[j*ntv+i]);
                        fprintf(DebEnv.t_logfile[mythid],"\n");
                     }
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  // Compute w = V^T * g, store w temporarily in work
                  cblas_dgemv(CblasColMajor,CblasNoTrans,rank_BB,ntv,1.0,VT,
                              ntv,&(g_scr[ind_g]),1,0.0,work,1);

                  // Scale w with the inverse of singular values
                  for (int i = 0; i < rank_BB; i++) work[i] /= SIGMA[i];

                  // Compute coef_P0 += V * w
                  cblas_dgemv(CblasColMajor,CblasNoTrans,nr_BB_loc,rank_BB,1.0,
                              &(BB_scr[ind_BB]),nr_BB_loc,work,1,1.0,
                              &(coef_P0[istart_patt]),1);
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  if (DEBUG){
                     fprintf(DebEnv.t_logfile[mythid],"\nDP:\n");
                     for (int j = 0; j < nr_BB_loc; j++)
                       fprintf(DebEnv.t_logfile[mythid]," %15.6e",coef_P0[istart_patt+j]);
                     fprintf(DebEnv.t_logfile[mythid],"\n");
                     fflush(DebEnv.t_logfile[mythid]);
                  }
                  //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                  if (true) { //@@@@ DA VALUTARE SE METTERLO TRA I PARAMETRI
                     // Pad BB_scr with zeroes
                     int npad = (ntv-rank_BB)*nr_BB_loc;
                     for (int k = 0; k < npad; k++) BB_scr[ind_BB+nr_BB_loc*rank_BB+k] = 0.0;
                  }

                  // Compute the number of entries of Q
                  nnz_BB = nr_BB_loc*ntv;

               }

               // Update pointers in BB, g and tau
               ind_BB += nnz_BB;
               ind_g += ntv;
               ind_tau += ntv;

            }

         }

      } // End loop over columns
      exit_loop_icol: ;

      // Free private workspace
      free(work);
      free(tau);
      free(VT);
      free(SIGMA);

      // Exit point
      exit_pragma: ;

   } // End of parallel region
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   //of = fopen("P0_dopo","w");
   //for (int i = 0; i < iat_patt[nn]; i++) fprintf(of,"%20.11e\n",coef_P0[i]);
   //fclose(of);
   //cout << "FATTO" << endl;
   //exit(0);
   //fclose(bbf);
   //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   // Free shared scratches
   free(c2glo);
   free(vec_g);

   return ierr;

}
