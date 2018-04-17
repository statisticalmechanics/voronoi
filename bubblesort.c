#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

//void BubbleSort(double *rx[],double *ry[],double *rz[],double *rij2[],int *tag[],int ncan)
void BubbleSort(int ncan)
{
 int i,itop,i1,tagi;
 double rxi,ryi,rzi,rij2i;

 int change;

 itop = ncan-1;

do
{
 change = FALSE;
 for(i=0;i<itop;i++)
 {
  i1 = i+1;

  if(dRij2[i] > dRij2[i1])
  {
   rxi =  dRX[i];
   ryi =  dRY[i];
   rzi =  dRZ[i];
   rij2i = dRij2[i]; 
   tagi = Tag[i];

   dRX[i] = dRX[i1];
   dRY[i] = dRY[i1];
   dRZ[i] = dRZ[i1];
   dRij2[i] = dRij2[i1]; 
   Tag[i] = Tag[i1];

   dRX[i1] = rxi;
   dRY[i1] = ryi;
   dRZ[i1] = rzi;
   dRij2[i1] = rij2i; 
   Tag[i1] = tagi;

   change = TRUE;
  }
 }//end loop i
 itop--;
}
while(change==TRUE && itop > 0);
//while(change==TRUE);
  
 return;
}
