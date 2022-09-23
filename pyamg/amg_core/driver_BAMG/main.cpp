#include <iostream>  // to use: cout,endl
#include <stdlib.h>  // to use: exit
#include <fstream>   // to use: ifstream,ofstream
#include <libgen.h>  // to use: basename
#include "omp.h"     // to use: omp_get_wtime
using namespace std;

#include "readMat.h"
#include "readMIS.h"
#include "prol_parm.h"
#include "cptBAMGProl.h"
#define charStrLen 2000

// MAIN PROGRAM
int main(int argc, const char* argv[]){

FILE *paramFile;

// Main input parameters
int stronly;
int np,itmax_vol,dist_min,dist_max,mmax;
double toll_vol,condfact,mincond,eps;
char SoC_file[charStrLen]="", MIS_file[charStrLen]="", TV_file[charStrLen]="";

// Other variables
int ierr, kk;
int nn_S, nt_S;
int *iat_S, *ja_S;
double *coef_S_tmp;
int *coef_S;
int nn_C, nn_F;
int *fcnodes;
int ntv;
double **TV;
double *TVbuf;
prol_parm params;
int *vecstart;
int nn_I, nc_I, nt_I;
int *iat_I, *ja_I;
double *coef_I;

char line[1024];

   // Check arguments
   if (argc < 2) {
      printf("Too few arguments.\n Usage: [param file] \n");
      exit(1);
   }

   // Read main inputs
   if ( (paramFile = fopen(argv[1], "r")) != NULL) {
      fgets(line, sizeof line, paramFile); sscanf(line, "%d",  &np);
      fgets(line, sizeof line, paramFile); sscanf(line, "%d",  &stronly);
      fgets(line, sizeof line, paramFile); sscanf(line, "%lf", &condfact);
      fgets(line, sizeof line, paramFile); sscanf(line, "%lf", &mincond);
      fgets(line, sizeof line, paramFile); sscanf(line, "%d",  &itmax_vol);
      fgets(line, sizeof line, paramFile); sscanf(line, "%lf", &toll_vol);
      fgets(line, sizeof line, paramFile); sscanf(line, "%lf", &eps);
      fgets(line, sizeof line, paramFile); sscanf(line, "%d",  &dist_min);
      fgets(line, sizeof line, paramFile); sscanf(line, "%d",  &dist_max);
      fgets(line, sizeof line, paramFile); sscanf(line, "%d",  &mmax);
      fgets(line, sizeof line, paramFile); sscanf(line, "%s",  SoC_file);
      fgets(line, sizeof line, paramFile); sscanf(line, "%s",  MIS_file);
      fgets(line, sizeof line, paramFile); sscanf(line, "%s",  TV_file);
      fclose(paramFile);
   } else {
      printf(" Error while opening file %s\n", argv[1]);
      exit(0);
   }

   // Dump input
   cout << "np :           " << np << endl;
   cout << "stronly:       " << stronly << endl;
   cout << "condfact:      " << condfact << endl;
   cout << "mincond:       " << mincond << endl;
   cout << "itmax_vol:     " << itmax_vol << endl;
   cout << "toll_vol:      " << toll_vol << endl;
   cout << "eps:           " << eps << endl;
   cout << "dist_min:      " << dist_min << endl;
   cout << "dist_max:      " << dist_max << endl;
   cout << "mmax:          " << mmax << endl;
   cout << "SoC_file:      " << SoC_file << endl;
   cout << "MIS_file:      " << MIS_file << endl;
   cout << "TV_file:       " << TV_file << endl;
   cout << endl;

   // Store input in the structure
   params.stronly = stronly;
   params.condfact = condfact;
   params.mincond = mincond;
   params.itmax_vol = itmax_vol;
   params.toll_vol = toll_vol;
   params.eps = eps;
   params.dist_min = dist_min;
   params.dist_max = dist_max;
   params.mmax = mmax;

   // Read the SoC matrix
   ierr = readCSRmat(&nn_S, &kk, &nt_S, &iat_S, &ja_S, &coef_S_tmp, SoC_file, 0);
   if (ierr != 0) exit(ierr);
   coef_S = new int [nt_S]();
   for (int i = 0; i < nt_S; i++) coef_S[i] = 1;

   cout << "SoC matrix read" << endl;
   /////////////////////////////
   //ofstream xfile("XXX");
   //for (int i=0; i<nt_S; i++) xfile << coef_S[i] << "  " << coef_S_tmp[i] << "  " << (coef_S_tmp[i] >0.0) << endl;
   //xfile.close();
   /////////////////////////////

   // Read the F/C partition
   ierr = readMIS(nn_S, nn_F, nn_C, fcnodes, MIS_file);
   if (ierr != 0) exit(ierr);
   cout << "F/C partition read" << endl;
   cout << "coarse nodes " << nn_C << endl;

   // Read the test space
   // Open the input file
   ifstream fileTV(TV_file);
   // Read header
   int k = 0;
   std::string str;
   fileTV >> str >> str >> ntv >> str;
   //if (kk != nn_S) exit(1);
   // Allocate and read TV matrix
   TV   = (double**) malloc( nn_S*sizeof(double*) );
   TVbuf = (double*) malloc( (ntv*nn_S)*sizeof(double) );
   for (int i = 0; i < nn_S; i++){
      TV[i] = &(TVbuf[k]);
      k += ntv;
      for (int j = 0; j < ntv; j++) fileTV >> TV[i][j];
   }
   // Close the input fileTV
   fileTV.close();
   cout << "TSPACE read " << ntv << endl;

   // Create data structure to manage omp
   vecstart = new int [np+1]();
   {  int bsize = nn_S/np;
      int resto = nn_S%np;
      vecstart[0] = 0;
      for (int i = 1; i <= resto; i++){
         vecstart[i] = vecstart[i-1] + bsize + 1;
      }
      for (int i = resto+1; i <= np; i++){
         vecstart[i] = vecstart[i-1] + bsize;
      }
   }

   for (int i = 0; i <= np; i++){ cout << vecstart[i] << " ";}

   // Get wall time
   double prol_time = omp_get_wtime();
   int * c_mark = (int*) malloc( nn_S*sizeof(int) );

   // Call the prolongation construction
   int level = 1;
   int verbosity = 1;
   double maxcond = 1.e13;
   double maxrownrm = 5.0;
   nn_I = nn_S;
   nc_I = nn_C;
   nt_I = nn_I*mmax;
   iat_I = (int*) malloc( (nn_I+1)*sizeof(int) );
   ja_I = (int*) malloc( nt_I*sizeof(int) );
   coef_I = (double*) malloc( nt_I*sizeof(double) );

   ierr = cptBAMGProl( level,
                       verbosity,
                       itmax_vol,
                       dist_min,
                       dist_max,
                       mmax,
                       maxcond,
                       maxrownrm,
                       toll_vol,
                       eps,
                       nn_S,
                       iat_S, nn_S+1,
                       ja_S, nt_S,
                       ntv,
                       fcnodes, nn_S,
                       TVbuf, ntv*nn_S,
                       nn_I,
                       nt_I,
                       iat_I, nn_I+1,
                       ja_I, nt_I,
                       coef_I, nt_I,
                       c_mark, nn_I );

   if (ierr != 0) exit(ierr);

   // Print the time required for prolongation
   prol_time = omp_get_wtime() - prol_time;
   cout << endl;
   cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   cout << "Time for building Prolongation [s]: " << fixed << prol_time << endl;
   cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
   cout << endl;

   // Print prolongation on output file
   cout << "Printing Prolongation Matrix" << endl;
   FILE *ofile = fopen("mat_P.csr","w"); if (!ofile) exit(1);
   for (int i = 0; i < nn_I; i++){
       for (int j = iat_I[i]; j < iat_I[i+1]; j++){
          fprintf(ofile,"%10d %10d %25.15e\n",i+1,ja_I[j]+1,coef_I[j]);
       }
   }
   fclose(ofile);

   // Delete scratch vectors
   delete [] coef_S;
   delete [] TV;
   delete [] TVbuf;
   delete [] vecstart;

   exit(0);
}
