/*******************************************
 * 3D Voronoi tessellation
 * binary/ternary mixtures
 * radical plane method
 * insize cubic box under periodic boundary condition
 * input particle coorinates and sizes
 * output neigbor list, vertices,edges
 * modified from Allen and Tildesley, 1987
 * Kai Zhang, Yale University, 2013
 *******************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
//#include <time.h>

int main(int argc, char *argv[])
{
 char filename[20];
 FILE *fpmovie;
 FILE *fp;
 FILE *fpvoro;
 FILE *fpvertex;
 FILE *fpcluster;
 FILE *fpedge;
 int i,j,k;
 int iv,jv;//vertex index
 double xij,yij,zij,rij2,sigmaij;
 int nbond,ibond[100000],jbond[100000];

 int can,Ncan; // candidates index and number

 int line,snapshot;
 char AtomID;
 double T;

 sscanf(argv[1],"%d",&JobIndex);

 ReadInput();
 
 sprintf(filename,"movie_%d.xyz",JobIndex);
 fpmovie = fopen(filename,"r");
 sprintf(filename,"bond_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 
for(snapshot=0;snapshot<Nsample;snapshot++)
{
  printf("snapshot %d\n",snapshot);
 for(line=0;line<=N+1;line++)
 {
  if(line>1)
  fscanf(fpmovie,"%c %lf %lf %lf\n",&AtomID,&RX[line-2],&RY[line-2],&RZ[line-2]);
  else if(line == 0)
  fscanf(fpmovie,"%d\n",&N);
  else
  fscanf(fpmovie,"%lf %lf\n",&L,&T);
 }
 V=CUBIC(L);

 for(i=0;i<N;i++)
 {
  NNeighbor[i] = 0;
  for(j=0;j<MAX_Neighbors;j++)
   NeighborList[i][j] = -1;
 }

 
 for(i=0;i<N;i++)
 {
  can = 0;
  for(j=0;j<N;j++)
  {
   if(j != i)
   {
    xij = RX[j]-RX[i]; // Rij from center i to neighbor j
    yij = RY[j]-RY[i];
    zij = RZ[j]-RZ[i];
    xij = xij - L*round(xij/L);
    yij = yij - L*round(yij/L);
    zij = zij - L*round(zij/L);
    rij2 = SQR(xij) + SQR(yij) + SQR(zij);
    sigmaij = (sigma[i]+sigma[j])/2.;

    if(rij2 < rc2*SQR(sigmaij))
    {
     dRX[can] = xij;
     dRY[can] = yij;
     dRZ[can] = zij;
     dRij2[can] = rij2;
     Tag[can] = j;
 //   if(i==0) printf("(%lf %lf %lf)\n",RX[j],RY[Tag[can]],RZ[Tag[can]]);
 //   if(i==0) printf("can = %d Tag = %d rij = %lf\n",can,Tag[can],rij2);
     can++;
    }//end if rij < rcutoff
   }//end if j!= i
  }//end loop j
  Ncan = can; // number of candidates of center i

 // printf("%d candidate neighbors %d\n",i,Ncan);
  //BubbleSort(&(dRX)[],&(dRY[]),&(dRZ[]),&(dRij2[]),&(Tag[]),Ncan); // sort candidates by distance dRij2
  BubbleSort(Ncan); // sort candidates by distance dRij2

  Work(i,Ncan);

  for(can=0;can<Ncan;can++)
  {
   if(Edges[can] > 0)
   {
    NeighborList[i][NNeighbor[i]] = Tag[can];
    NNeighbor[i]++ ;

   // if(Tag[can]>i) printf("%d %d\n",i,Tag[can]);
   }
  }

   sprintf(filename,"tessellation_%d.dat",JobIndex);
   fpvoro=fopen(filename,"a+");
   fprintf(fpvoro,"i = %d NF = %d NV = %d NE = %d\n",i,NFaces,NVertices,NEdges);
   fclose(fpvoro);

 // if(SoluteCheck==TRUE)
  {
   sprintf(filename,"cluster_%d.dat",JobIndex);
   fpcluster=fopen(filename,"a+");
   fprintf(fpcluster,"i = %d NF3 = %d NF4 = %d NF5 = %d NF6 = %d NF7 = %d NF8 = %d NF9 = %d check = %d\n",i,NF3,NF4,NF5,NF6,NF7,NF8,NF9,NF3+NF4+NF5+NF6+NF7+NF8+NF9 == NFaces);
   fclose(fpcluster);
  }
  
  sprintf(filename,"vertex_%d.xyz",JobIndex);
  fpvertex=fopen(filename,"w");
  for(iv=0;iv<NVertices;iv++)
   fprintf(fpvertex,"H %lf %lf %lf\n",RX[i]+dRVX[iv],RY[i]+dRVY[iv],RZ[i]+dRVZ[iv]);
  fclose(fpvertex);
  

   
  sprintf(filename,"edge_%d.dat",JobIndex);
  if(JobIndex == 0) 
  fpedge=fopen(filename,"a+"); // make "w" to print only current; "a+" to print all
  else
  fpedge=fopen(filename,"w"); // make "w" to print only current; "a+" to print all

  fprintf(fpedge,"draw color red\n");
  for(iv=0;iv<NVertices-1;iv++)
  for(jv=iv+1;jv<NVertices;jv++)
  {
   if(ShareEdge(iv,jv)==TRUE) 
   fprintf(fpedge,"draw line \"%lf %lf %lf\" \"%lf %lf %lf\" width 2\n",RX[i]+dRVX[iv],RY[i]+dRVY[iv],RZ[i]+dRVZ[iv],RX[i]+dRVX[jv],RY[i]+dRVY[jv],RZ[i]+dRVZ[jv]);
  }
  fclose(fpedge);

 }//end loop i

 nbond=0;
 for(i=0;i<N-1;i++)
 for(j=i+1;j<N;j++)
 {
  for(can=0;can<NNeighbor[i];can++)
  if(NeighborList[i][can]==j)
  {
   ibond[nbond]=i;
   jbond[nbond]=j;
   nbond++;
  }
 }
 
 fprintf(fp,"%d\n",nbond);
 for(i=0;i<nbond;i++)
   fprintf(fp,"%d %d\n",ibond[i],jbond[i]);
}//end loop over snapshots
 fclose(fpmovie);
 fclose(fp);

 /*
 sprintf(filename,"movie_%d.xyz",JobIndex);
 fp=fopen(filename,"w");
 Writemovie(fp);
 fclose(fp);
 */

 printf("*******************************************************************************\n");
 return 0;
}
