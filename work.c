#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void Work(int I0,int ncan)
{
 int i,j,k,l;
 int nver; // vertex index
 double tol;

 int OK; 
 double ai,bi,ci,di; //ax+by+cz = d
 double aj,bj,cj,dj;
 double ak,bk,ck,dk;
 double det;

 double ab,bc,ca,da,db,dc;
 double vxijk,vyijk,vzijk; // vertex coordinate relative to center particle

 tol = 1.E-6;

 if(ncan < 4)
 {
  printf("less than 4 neighbors %d\n",ncan);
 // exit(1);
  return;
 }

 nver=0;
  


 for(i=0;i<ncan-2;i++)
 {
 ai = dRX[i];
 bi = dRY[i];
 ci = dRZ[i];
 di = 0.5*(dRij2[i]+(SQR(sigma[I0]/2.)-SQR(sigma[Tag[i]]/2.))); // differ by -1* from Allen and Tildesley
 for(j=i+1;j<ncan-1;j++)
 {
 aj = dRX[j];
 bj = dRY[j];
 cj = dRZ[j];
 dj = 0.5*(dRij2[j]+(SQR(sigma[I0]/2.)-SQR(sigma[Tag[j]]/2.)));

 ab = ai*bj-aj*bi;
 bc = bi*cj-bj*ci;
 ca = ci*aj-cj*ai;
 da = di*aj-dj*ai;
 db = di*bj-dj*bi;
 dc = di*cj-dj*ci;

 for(k=j+1;k<ncan;k++)
 {
 ak = dRX[k];
 bk = dRY[k];
 ck = dRZ[k];
 dk = 0.5*(dRij2[k]+(SQR(sigma[I0]/2.)-SQR(sigma[Tag[k]]/2.)));

 det = ak*bc+bk*ca+ck*ab;

 if(fabs(det) > tol)
 {
 vxijk = (dk*bc-bk*dc+ck*db)/det;
 vyijk = (ak*dc+dk*ca-ck*da)/det;
 vzijk = (-ak*db+bk*da+dk*ab)/det;

 OK = TRUE;
 l = 0;
 do
 {
  if(l!=i && l!=j && l!=k)
  {
   if(dRX[l]*vxijk+dRY[l]*vyijk+dRZ[l]*vzijk <= 0.5*(dRij2[l]+(SQR(sigma[I0]/2.)-SQR(sigma[Tag[l]]/2.)))) // if plane l is beyond vertex
   OK = TRUE;
   else 
   OK = FALSE;
  }
  l++;
 }
 while(OK==TRUE && l<ncan);


 if(OK==TRUE)
 {
  IVer[nver]=i;
  JVer[nver]=j;
  KVer[nver]=k;

  dRVX[nver] = vxijk;
  dRVY[nver] = vyijk;
  dRVZ[nver] = vzijk;

  nver++;
 }

 }// end if det != 0

 //printf("(%d %d %d)check nver=%d\n",i,j,k,nver);

 }// end loop k
 }// end loop j
 }// end loop i
  

 NVertices = nver;

 if(NVertices < 4)
 {
 printf("%d less than 4 vertices\n",NVertices);
 //exit(1);
 return;
 }

 for(i=0;i<ncan;i++)
  Edges[i] = 0;
 
 for(nver=0;nver<NVertices;nver++)
 {
  Edges[IVer[nver]]++;
  Edges[JVer[nver]]++;
  Edges[KVer[nver]]++;
 }

 NFaces=0;
  NF3=0;
  NF4=0;
  NF5=0;
  NF6=0;
  NF7=0;
  NF8=0;
  NF9=0;
 NEdges=0;

 SoluteCheck=TRUE;
 for(i=0;i<ncan;i++)
 {
  if(Edges[i]>0) 
  {
  NFaces++;
  if(identity[Tag[i]] == identity[I0])
  SoluteCheck=FALSE;
  }

  if(Edges[i]==3) NF3++;
  if(Edges[i]==4) NF4++;
  if(Edges[i]==5) NF5++;
  if(Edges[i]==6) NF6++;
  if(Edges[i]==7) NF7++;
  if(Edges[i]==8) NF8++;
  if(Edges[i]==9) NF9++;
  
  NEdges += Edges[i];
 }
  
 if(NEdges%2 !=0)
 {
  printf("noninteger number of edges 2NE=%d\n",NEdges);
 // exit(1);
 }
 
 NEdges = NEdges/2;

 if(NVertices-NEdges+NFaces != 2)
 {
  printf("NV=%d NE=%d NF=%d Euler Error: degeneracy?\n",NVertices,NEdges,NFaces);
 }
 
 if(NF3+NF4+NF5+NF6+NF7+NF8+NF9 != NFaces)
  printf("NF3+NF4+NF5+NF6+NF7+NF8+NF9 =%d != %d = NF\n",NF3+NF4+NF5+NF6+NF7+NF8+NF9,NFaces);




 return;
}


int ShareEdge(int i,int j) // whether two vertices connected edge
{
 int check;

 int ai,bi,ci; // ai<bi<ci
 int aj,bj,cj;

 check=FALSE;// don't share

 ai=MIN(MIN(IVer[i],JVer[i]),KVer[i]);
 ci=MAX(MAX(IVer[i],JVer[i]),KVer[i]);
 bi = IVer[i]+JVer[i]+KVer[i]-ai-ci; 
 
 aj=MIN(MIN(IVer[j],JVer[j]),KVer[j]);
 cj=MAX(MAX(IVer[j],JVer[j]),KVer[j]);
 bj = IVer[j]+JVer[j]+KVer[j]-aj-cj; 

 if((ai==aj && bi==bj) || (ai==bj && bi==cj) || (ai==aj && bi==cj) \
     || (bi==bj && ci==cj) || (bi==aj && ci==bj) || (bi==aj && ci==cj) \
     || (ai==aj && ci==cj) || (ai==aj && ci==bj) || (ai==bj && ci==cj))
 {
 check=TRUE;
 }

// printf("(%d %d %d) and (%d %d %d) check=%d\n",ai,bi,ci,aj,bj,cj,check);
 return check;
}


