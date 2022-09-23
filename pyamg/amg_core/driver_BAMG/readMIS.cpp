#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <fstream>  // to use: ifstream
/////////////////////
#include <iostream>
using namespace std;
/////////////////////

int readMIS(int const nn, int & nn_F, int & nn_C, int *& fcnodes, char * fname)
{
   // Open MIS file
   printf("MIS file %s\n",fname);
   ifstream file(fname);
   // Skip two lines
   fcnodes = (int*) malloc(nn * sizeof(int));
   nn_C = 0;
   for (int i = 0; i < nn; i++){
      file >> fcnodes[i];
      if (fcnodes[i] >= 0) {
         nn_C++;
      }
   }
   nn_F = nn - nn_C;

   // close the input file
   file.close();

   return 0;
}
