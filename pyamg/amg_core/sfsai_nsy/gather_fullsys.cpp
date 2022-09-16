void gather_fullsys(const int mrow, const int *const vecinc, const int nn,
                    const int *const iat, const int *const ja, const double *const coef,
                    double &diag_entry, double *full_A, double *rhs_L, double *rhs_U,
                    bool &null_L, bool &null_U){

   // Loop over the rows(+1) of full_A
   null_L = true;
   null_U = true;
   int ind_row = 0;
   for (int i = 0; i < mrow; i++){
      // Load i-th row of full_A exploring the irow-th row of matrix A
      int ii = 0;
      int irow = vecinc[i];
      int jj = iat[irow];
      int endrow = iat[irow+1];
      while (ii < mrow){
         // Check ja[jj] >= vecinc[ii]
         while (ja[jj] < vecinc[ii]){
            jj++;
            if (jj == endrow){
               // End of row reached, set to zero the remaining entries and go to next row
               for (int k = ii; k < mrow; k++) full_A[ind_row+k] = 0.0;
               rhs_U[i] = 0.0;
               goto next_row;
            }
         }
         if (vecinc[ii] == ja[jj]){
            // Add this entry to full_A
            full_A[ind_row+ii] = coef[jj];
            ii++;
         } else {
            // This entry is null
            full_A[ind_row+ii] = 0.0;
            ii++;
         }
      } // End external while loop
      // Set the i-th entry of rhs_U
      while(ja[jj] < vecinc[ii]){
         jj++;
         if (jj == endrow){
            rhs_U[i] = 0.0;
            goto next_row;
         }
      }
      if (vecinc[ii] == ja[jj]){
         // Add this entry to rhs_U
         rhs_U[i] = -coef[jj];
         null_U = false;
      } else {
         // This entry is null
         rhs_U[i] = 0.0;
      }
      next_row:;
      ind_row += mrow;
   } // End Row loop

   // Load rhs_L only exploring the vencinc[mrow]-th row of matrix A
   int ii = 0;
   int irow = vecinc[mrow];
   int jj = iat[irow];
   int endrow = iat[irow+1];
   while (ii < mrow){
      // Check ja[jj] >= vecinc[ii]
      while (ja[jj] < vecinc[ii]){
         jj++;
         if (jj == endrow){
            // End of row reached, set to zero the remaining entries and go to next row
            for (int k = ii; k < mrow; k++) rhs_L[k] = 0.0;
            goto end_rhs_L;
         }
      }
      if (vecinc[ii] == ja[jj]){
         // Add this entry to rhs_L
         rhs_L[ii] = -coef[jj];
         null_L = false;
         ii++;
      } else {
         // This entry is null
         rhs_L[ii] = 0.0;
         ii++;
      }
      end_rhs_L:;
   } // End external while loop

   // Retrieve diagonal entry of A
   while (jj < endrow){
      if (ja[jj] == irow){
         diag_entry = coef[jj];
         return;
      }
      jj++;
   }
   diag_entry = 0.0;

}
