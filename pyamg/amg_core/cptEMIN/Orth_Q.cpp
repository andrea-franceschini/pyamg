#include "omp.h"
#include "cblas.h"

// Performs the orthogonalization v_out = v_in - Q*QT*v_in
void Orth_Q(const int np, const int nn, const int nn_K, const int ntv,
            const int* __restrict__ iat_patt, const int* __restrict__ ja_patt,
            const double* __restrict__ mat_Q, const double* __restrict__ v_in,
            double* __restrict__ v_ntv, double* __restrict__ v_out){

   #pragma omp parallel num_threads(np)
   {
      // Select the local scratch for v_ntv
      int myid = omp_get_thread_num();
      double *my_v_ntv = &(v_ntv[myid*ntv]);
      #pragma omp for
      for (int icol = 0; icol < nn; icol++){
         int istart = iat_patt[icol];
         int iend = iat_patt[icol+1];
         int nr_loc = iend - istart;
         if (nr_loc > 0){
            int ind_Q = istart*ntv;

            // Perform product by QT: my_v_ntv <-- QT*v_in;
            cblas_dgemv(CblasColMajor,CblasTrans,nr_loc,ntv,1.0,&(mat_Q[ind_Q]),nr_loc,
                        &(v_in[istart]),1,0.0,my_v_ntv,1);

            // Perform product by Q: v_out <-- Q*my_v_ntv;
            cblas_dgemv(CblasColMajor,CblasNoTrans,nr_loc,ntv,1.0,&(mat_Q[ind_Q]),nr_loc,
                        my_v_ntv,1,0.0,&(v_out[istart]),1);

            // Compute v_out <-- v_in - v_out
            for (int j = istart; j < iend; j++) v_out[j] = v_in[j] - v_out[j];
         }
      }
   }

}
