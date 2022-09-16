#pragma once

int cpt_nsy_sfsai_coef(const double tau, const int mmax, const int nn,
                       const int *const iat_A, const int *const ja_A,
                       const double *const coef_A, const int *const iat_F,
                       const int *const ja_F,  int &nt_FL, int &nt_FU,
                       double *coef_FL, double *coef_FU);
