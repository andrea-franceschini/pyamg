typedef int iReg;
typedef int iExt;
typedef double rExt;

void UinvDP_spmat(const iReg nthreads, const iReg nrows_A, const iExt* iat_A,
                  const iReg *ja_A, rExt *coef_A, const iReg ncols_B, const iExt *iat_B,
                  const iReg *ja_B, const rExt *coef_B, iReg* WNALL);
