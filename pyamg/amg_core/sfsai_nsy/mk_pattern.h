#pragma once

int mk_pattern(const int kpow, const double tau_pref, const int nnzr_max, const int nn_A,
               const int nt_A, const int *const iat_A, const int *const ja_A,
               const double *const coef_A, int &mmax, int * iat_Pout, const int ja_Pout_size,
               int * ja_Pout);
