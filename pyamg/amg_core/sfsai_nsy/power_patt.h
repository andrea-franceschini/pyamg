#pragma once

int power_patt(const int kpow, const int nnzr_max, const int nn, const int *const iat,
               const int *const ja, int &mmax, int * iat_Pout, const int ja_Pout_size,
               int * ja_Pout);
