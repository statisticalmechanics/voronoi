#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

void Writemovie(FILE *FilePtr)
{
  int i;
  
  fprintf(FilePtr,"%d\n",N);
  fprintf(FilePtr,"%lf\n",L);

  for(i=0;i<N;i++)
  {
     if(identity[i] == 1)	
      fprintf(FilePtr,"%s\t","N");
     else if(identity[i] == 2)	
      fprintf(FilePtr,"%s\t","S");
     else
      fprintf(FilePtr,"%s\t","O");

     fprintf(FilePtr,"%lf\t",RX[i]);
     fprintf(FilePtr,"%lf\t",RY[i]);
     fprintf(FilePtr,"%lf\n",RZ[i]);
  }
}
