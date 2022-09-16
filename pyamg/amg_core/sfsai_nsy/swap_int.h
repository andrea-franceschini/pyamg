#pragma once

inline void swap_int(int &i1, int &i2){
   int tmp = i1;
   i1 = i2;
   i2 = tmp;
   return;
}
