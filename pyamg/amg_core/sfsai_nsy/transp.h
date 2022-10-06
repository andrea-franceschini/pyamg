int transp_patt(const int np, const int nrows, const int ncols,
                const int* __restrict__ iat, const int* __restrict__ ja,
                int*& __restrict__ iat_T, int*& __restrict__ ja_T);

int transp_csrmat(const int np, const int nrows, const int ncols,
                  const int* __restrict__ iat, const int* __restrict__ ja,
                  const double* __restrict__ coef, int*& __restrict__ iat_T,
                  int*& __restrict__ ja_T, double*& __restrict__ coef_T);
