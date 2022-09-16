
//----------------------------------------------------------------------------------------

/*
 * Sorts an integer array x1 in such a way that x1(i) <= x1(i+1)
 *
 * > @brief Compute the gradient of the Kaporin Number
 *
 * > @param[in]     n  # of components of x1
 * ! @param[inout]  x1 array of integers
 *
 * > @author Carlo Janna
 * ! @author Giovanni Isotton
 *
 * > @version 1.0
 *
 * > @date april 2019
 *
 * > @par License
 * ! This program is intended for private use only and can not be distributed
 * ! elsewhere without authors' consent
 *
 */

#include "swap_int.h"

void iheapsort(int *x1, int n){

for (int node = 2; node < n+1; node ++){
   int i = node;
   int j = i/2;
   while( x1[j-1] < x1[i-1] ){
      swap_int(x1[j-1],x1[i-1]);
      i = j;
      j = i/2;
      if (i == 1) break;
   }
}

for (int i = n; i > 1; i --){
   swap_int(x1[i-1],x1[0]);
   int k = i - 1;
   int ik = 1;
   int jk = 2;
   if (k >= 3){
      if (x1[2] > x1[1]) jk = 3;
   }
   bool cont_cycle = false;
   if (jk <= k){
      if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
   }
   while (cont_cycle){
      swap_int(x1[jk-1],x1[ik-1]);
      ik = jk;
      jk = ik*2;
      if (jk+1 <= k){
         if (x1[jk] > x1[jk-1]) jk = jk+1;
      }
      cont_cycle = false;
      if (jk <= k){
         if (x1[jk-1] > x1[ik-1]) cont_cycle = true;
      }
   }
}

}

//----------------------------------------------------------------------------------------
