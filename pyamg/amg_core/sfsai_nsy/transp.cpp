#include "omp.h"
#include <stdlib.h>
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
//@@@@@@@@@@@@@@@@@@@@@@@@@@
//#include <iostream>
//using std::endl;
//using std::cout;
//@@@@@@@@@@@@@@@@@@@@@@@@@@

//----------------------------------------------------------------------------------------

// Counts the number of non-zeroes assigned to a thread in each row
void count_rowterms( const int nequ, const int nterm, const int* __restrict__ ja,
                     int* __restrict__ WI ){

   // Initialize WI
   for ( int i = 0; i < nequ; i++ ) {
      WI[i] = 0;
   }

   // Count non-zeroes
   for ( int i = 0; i < nterm; i++ ) {
      WI[ja[i]]++;
   }

}

//----------------------------------------------------------------------------------------

// Computes the global iat_T and update WI pointers.
void mkiat_Tglo (const int myid, const int nrows, const int nequ, const int nthreads,
                 const int firstrow, int** __restrict__ WI,
                 const int* __restrict__  nnz, int* __restrict__ iat_T  ){

   // Compute the number of non-zeros belonging to previous threads
   int ntprec = 0;
   for ( int i = 0; i < myid; i++ ){
      ntprec += nnz[i];
   }

   // Update iat_T
   for ( int i = 0; i < nrows; i++ ){
      iat_T[i] += ntprec;
   }

   // Account for (nrows+1)-th component for last threads
   if ( myid + 1 == nthreads ) {
      iat_T[nrows] = ntprec + nnz[myid];
   }

   // Update WI1
   for ( int i = firstrow; i < firstrow+nrows; i++ ){
      WI[0][i] = iat_T[i-firstrow];
   }
   for ( int j = 1; j < nthreads; j++ ){
      for ( int i = firstrow; i < firstrow+nrows; i++ ){
         WI[j][i] +=  WI[j-1][i];
      }
   }

}

//----------------------------------------------------------------------------------------

// Each thread computes the pointers to the beginning of each row of its part of the matrix
void mkiat_Tloc (const int nrows, const int nequ, const int nthreads, const int firstrow,
                 int** __restrict__ WI, int* __restrict__ iat_T, int& nnz  ){

   // Perform reduction
   for ( int i = 0; i < nrows; i++ ) {
      iat_T[i] = WI[1][firstrow+i];
   }
   for ( int j = 2; j < nthreads+1; j++ ) {
      for ( int i = 0; i < nrows; i++ ) {
         iat_T[i] += WI[j][firstrow+i];
      }
   }

   // Compute iat_T
   int tmp = iat_T[0];
   iat_T[0] = 0;
   for ( int i = 1; i < nrows; i++ ) {
      int tmp2 = iat_T[i];
      iat_T[i] = iat_T[i-1] + tmp;
      tmp = tmp2;
   }

   // Compute the number of non-zeroes
   nnz = iat_T[nrows-1] + tmp;

}

//----------------------------------------------------------------------------------------

// Transpose matrix indices
void mvjcols(const int firstrow, const int nrows, const int nequ, const int nterm,
             const int nterm_T, const int* __restrict__ iat,
             const int* __restrict__ ja, int* __restrict__ ja_T,
             int* __restrict__ punt){

   // Transpose local stripe of F
   int shift = firstrow;
   for ( int i = 0; i < nrows; i++ ) {
      for ( int j = iat[i]; j < iat[i+1]; j++ ) {
         int irow = ja[j];
         ja_T[punt[irow]] = i + shift;
         punt[irow]++;
      }
   }

}

//----------------------------------------------------------------------------------------

// Transpose matrix indices and entries
void mvcoef(const int firstrow, const int nrows, const int nequ, const int nterm,
            const int nterm_T, const int* __restrict__ iat,
            const int* __restrict__ ja, const double* __restrict__ coef,
            int* __restrict__ ja_T, double* __restrict__ coef_T, int* __restrict__ punt){

   // Transpose local stripe of F
   int shift = firstrow;
   for ( int i = 0; i < nrows; i++ ) {
      for ( int j = iat[i]; j < iat[i+1]; j++ ) {
         int irow = ja[j];
         ja_T[punt[irow]] = i + shift;
         coef_T[punt[irow]] = coef[j];
         punt[irow]++;
      }
   }

}

//----------------------------------------------------------------------------------------

// Transpose a csr pattern
int transp_patt(const int np, const int nrows, const int ncols,
                const int* __restrict__ iat, const int* __restrict__ ja,
                int*& __restrict__ iat_T, int*& __restrict__ ja_T){

   // Init error code
   int ierr = 0;

   // Set number of threads
   int nthreads = MIN(np,MIN(nrows,ncols));

   // Allocate the transposed pattern
   int nterm = iat[nrows];
   iat_T = (int*) malloc((ncols+1) * sizeof(int));
   ja_T  = (int*) malloc((nterm) * sizeof(int));
   if (iat_T == nullptr || ja_T == nullptr) return ierr = 1;

   // Allocate scratches
   int* WI1_buf = (int*) malloc(((nthreads+1)*nrows) * sizeof(int));
   int** WI1 = (int**) malloc((nthreads+1) * sizeof(int*));
   int* WI2 = (int*) malloc(nthreads * sizeof(int));
   int* vecstart   = (int*) malloc((nthreads+1) * sizeof(int));
   int* vecstart_T = (int*) malloc((nthreads+1) * sizeof(int));
   if (WI1_buf == nullptr || WI1 == nullptr || WI2 == nullptr ||
       vecstart == nullptr || vecstart_T == nullptr) return ierr = 2;
   for (int i = 0; i <= nthreads; i++) WI1[i] = &(WI1_buf[i*nrows]);

   // Set pointer to the first unknown (row) of the omp thread
   int resto   = nrows%nthreads;
   int blksize = nrows/nthreads + 1;
   vecstart[0] = 0;
   for ( int i = 0; i < resto; i++) {
     vecstart[i+1] = vecstart[i] + blksize;
   }
   blksize -= 1;
   for ( int i = resto; i < nthreads; i++) {
     vecstart[i+1] = vecstart[i] + blksize;
   }

   // Set pointer to the first unknown (column) of the omp thread
   resto   = ncols%nthreads;
   blksize = ncols/nthreads + 1;
   vecstart_T[0] = 0;
   for ( int i = 0; i < resto; i++) {
     vecstart_T[i+1] = vecstart_T[i] + blksize;
   }
   blksize -= 1;
   for ( int i = resto; i < nthreads; i++) {
     vecstart_T[i+1] = vecstart_T[i] + blksize;
   }

   // Counts the number of non-zeroes assigned to each processor in each MatT row
   #pragma omp parallel num_threads(nthreads)
   {
      // Get thread ID
      int mythid = int( omp_get_thread_num() );

      // Retrieve thread info (row distribution)
      bool activeth = false;
      int firstrow = 0;
      int lastrow  = 0;
      int nrowsth  = 0;
      if ( mythid < nthreads ) {
         activeth = true;
         firstrow = vecstart[mythid];
         lastrow  = vecstart[mythid+1] - 1;
         nrowsth  = lastrow - firstrow + 1;
      }

      // Count
      int istart1 = iat[firstrow];
      int nnzth = iat[lastrow+1] - istart1;
      if ( activeth == true ) {
         if (nnzth > 0) {
            count_rowterms(ncols,nnzth,&ja[istart1],WI1[mythid+1]);
         } else {
            for ( int i = 0; i < ncols; i++ ) {
               WI1[mythid+1][i] = 0;
            }
         }
      }
      #pragma omp barrier
      /////////////////////////////////////////////
      //#pragma omp single
      //{
      //   for ( int ith = 0; ith < nthreads; ith++){
      //      for ( int i = 0; i < ncols; i++ ) {
      //         fprintf(t_logfile[ith]," %6d\n",WI1[ith+1][i]);
      //         //if ((i+1)%10 == 0) fprintf(t_logfile[ith],"\n");
      //         fflush(t_logfile[ith]);
      //      }
      //      fprintf(t_logfile[ith],"\n");
      //      fflush(t_logfile[ith]);
      //   }
      //}
      /////////////////////////////////////////////

      // Retrieve thread info (column distribution)
      bool activeth_T = false;
      int firstrow_T = 0;
      int lastrow_T  = 0;
      int nrowsth_T  = 0;
      if ( mythid < nthreads ) {
         activeth_T = true;
         firstrow_T = vecstart_T[mythid];
         lastrow_T  = vecstart_T[mythid+1] - 1;
         nrowsth_T  = lastrow_T - firstrow_T + 1;
      }

      // Each processor compute its part of iat_T
      if(activeth_T==true) mkiat_Tloc(nrowsth_T,ncols,nthreads,firstrow_T,
                                      WI1,&iat_T[firstrow_T],WI2[mythid]);
      #pragma omp barrier

      // Update iat_T and WI1 pointers
      if(activeth_T==true) mkiat_Tglo(mythid,nrowsth_T,ncols,nthreads,firstrow_T,
                                      WI1,WI2,&iat_T[firstrow_T]);
      #pragma omp barrier

      // Transpose coefficients
      if(activeth_T==true) mvjcols(firstrow,nrowsth,ncols,nterm,nterm,&iat[firstrow],
                                   ja,ja_T,WI1[mythid]);

   } // end parallel region

   // Free scratches
   free(WI1_buf);
   free(WI1);
   free(WI2);
   free(vecstart);
   free(vecstart_T);
   
   return ierr;

}

//----------------------------------------------------------------------------------------

// Transpose a csr matrix
int transp_csrmat(const int np, const int nrows, const int ncols,
                  const int* __restrict__ iat, const int* __restrict__ ja,
                  const double* __restrict__ coef, int*& __restrict__ iat_T,
                  int*& __restrict__ ja_T, double*& __restrict__ coef_T){

   // Init error code
   int ierr = 0;

   // Set number of threads
   int nthreads = MIN(np,MIN(nrows,ncols));

   // Allocate the transposed pattern
   int nterm = iat[nrows];
   iat_T = (int*) malloc((ncols+1) * sizeof(int));
   ja_T  = (int*) malloc((nterm) * sizeof(int));
   coef_T  = (double*) malloc((nterm) * sizeof(double));
   if (iat_T == nullptr || ja_T == nullptr || coef_T == nullptr) return ierr = 1;

   // Allocate scratches
   int* WI1_buf = (int*) malloc(((nthreads+1)*nrows) * sizeof(int));
   int** WI1 = (int**) malloc((nthreads+1) * sizeof(int*));
   int* WI2 = (int*) malloc(nthreads * sizeof(int));
   int* vecstart   = (int*) malloc((nthreads+1) * sizeof(int));
   int* vecstart_T = (int*) malloc((nthreads+1) * sizeof(int));
   if (WI1_buf == nullptr || WI1 == nullptr || WI2 == nullptr ||
       vecstart == nullptr || vecstart_T == nullptr) return ierr = 2;
   for (int i = 0; i <= nthreads; i++) WI1[i] = &(WI1_buf[i*nrows]);

   // Set pointer to the first unknown (row) of the omp thread
   int resto   = nrows%nthreads;
   int blksize = nrows/nthreads + 1;
   vecstart[0] = 0;
   for ( int i = 0; i < resto; i++) {
     vecstart[i+1] = vecstart[i] + blksize;
   }
   blksize -= 1;
   for ( int i = resto; i < nthreads; i++) {
     vecstart[i+1] = vecstart[i] + blksize;
   }

   // Set pointer to the first unknown (column) of the omp thread
   resto   = ncols%nthreads;
   blksize = ncols/nthreads + 1;
   vecstart_T[0] = 0;
   for ( int i = 0; i < resto; i++) {
     vecstart_T[i+1] = vecstart_T[i] + blksize;
   }
   blksize -= 1;
   for ( int i = resto; i < nthreads; i++) {
     vecstart_T[i+1] = vecstart_T[i] + blksize;
   }

   // Counts the number of non-zeroes assigned to each processor in each MatT row
   #pragma omp parallel num_threads(nthreads)
   {
      // Get thread ID
      int mythid = int( omp_get_thread_num() );

      // Retrieve thread info (row distribution)
      bool activeth = false;
      int firstrow = 0;
      int lastrow  = 0;
      int nrowsth  = 0;
      if ( mythid < nthreads ) {
         activeth = true;
         firstrow = vecstart[mythid];
         lastrow  = vecstart[mythid+1] - 1;
         nrowsth  = lastrow - firstrow + 1;
      }

      // Count
      int istart1 = iat[firstrow];
      int nnzth = iat[lastrow+1] - istart1;
      if ( activeth == true ) {
         if (nnzth > 0) {
            count_rowterms(ncols,nnzth,&ja[istart1],WI1[mythid+1]);
         } else {
            for ( int i = 0; i < ncols; i++ ) {
               WI1[mythid+1][i] = 0;
            }
         }
      }
      #pragma omp barrier
      /////////////////////////////////////////////
      //#pragma omp single
      //{
      //   for ( int ith = 0; ith < nthreads; ith++){
      //      for ( int i = 0; i < ncols; i++ ) {
      //         fprintf(t_logfile[ith]," %6d\n",WI1[ith+1][i]);
      //         //if ((i+1)%10 == 0) fprintf(t_logfile[ith],"\n");
      //         fflush(t_logfile[ith]);
      //      }
      //      fprintf(t_logfile[ith],"\n");
      //      fflush(t_logfile[ith]);
      //   }
      //}
      /////////////////////////////////////////////

      // Retrieve thread info (column distribution)
      bool activeth_T = false;
      int firstrow_T = 0;
      int lastrow_T  = 0;
      int nrowsth_T  = 0;
      if ( mythid < nthreads ) {
         activeth_T = true;
         firstrow_T = vecstart_T[mythid];
         lastrow_T  = vecstart_T[mythid+1] - 1;
         nrowsth_T  = lastrow_T - firstrow_T + 1;
      }

      // Each processor compute its part of iat_T
      if(activeth_T==true) mkiat_Tloc(nrowsth_T,ncols,nthreads,firstrow_T,
                                      WI1,&iat_T[firstrow_T],WI2[mythid]);
      #pragma omp barrier

      // Update iat_T and WI1 pointers
      if(activeth_T==true) mkiat_Tglo(mythid,nrowsth_T,ncols,nthreads,firstrow_T,
                                      WI1,WI2,&iat_T[firstrow_T]);
      #pragma omp barrier

      // Transpose coefficients
      if(activeth_T==true) mvcoef(firstrow,nrowsth,ncols,nterm,nterm,&iat[firstrow],
                                   ja,coef,ja_T,coef_T,WI1[mythid]);

   } // end parallel region

   // Free scratches
   free(WI1_buf);
   free(WI1);
   free(WI2);
   free(vecstart);
   free(vecstart_T);
   
   return ierr;

}

//----------------------------------------------------------------------------------------
