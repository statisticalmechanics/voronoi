#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void ReadInput(void)
{
 char filename[20];
 FILE *fp;
 int i;

 fp = fopen("input","r");
 fscanf(fp,"%d %d %d %d",&Nsample,&Nstart,&Nend,&samplefrequency);
 fscanf(fp,"%d %d %d %d",&N,&NA,&NB,&NC);
 fscanf(fp,"%lf",&V);
 fscanf(fp,"%lf %lf %lf",&sigmaA,&sigmaB,&sigmaC);
 fscanf(fp,"%lf %lf %lf",&sigmaAB,&sigmaBC,&sigmaCA);
 fscanf(fp,"%lf",&rc);
 fclose(fp);

 L=pow(V,1./3.);
 rc2 = rc*rc;

 //sprintf(filename,"position_%d",JobIndex);

// printf("0 (%lf %lf %lf)\n",RX[0],RY[0],RZ[0]);
// printf("1 (%lf %lf %lf)\n",RX[1],RY[1],RZ[1]);
// printf("2 (%lf %lf %lf)\n",RX[2],RY[2],RZ[2]);


 for(i=0;i<N;i++)
 {
  if(i<NA)
  {
  sigma[i] = sigmaA;
  identity[i] = 1;
  }
  else if(i >= NA && i< NA+NB)
  {
  sigma[i] = sigmaB;
  identity[i] = 2;
  }
  else
  {
  sigma[i] = sigmaC;
  identity[i] = 3;
  }
 }

 printf("N = %d NA = %d NB = %d NC = %d\n",N,NA,NB,NC);
 //printf("V = %lf L = %lf\n",V,L);
 printf("sigmaA = %lf sigmaB = %lf sigmaC = %lf\n",sigmaA,sigmaB,sigmaC);
 printf("sigmaAB = %lf sigmaBC = %lf sigmaCA = %lf\n",sigmaAB,sigmaBC,sigmaCA);
 printf("rcutoff = %lf\n",rc);

 return;
}
