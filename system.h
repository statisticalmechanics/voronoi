/******************************************
 * headfile containing most global variables
 ******************************************/
#include <stdio.h>

#define MAX_N 10000
#define MAX_Neighbors 100
#define MAX_Vertices 100
#define TRUE 1
#define FALSE 0

#define SQR(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))

extern int Nsample,Nstart,Nend,samplefrequency;

extern int JobIndex;
extern int N,NA,NB,NC; // number of particles in ternary mixtures: A,B,C. NA+NB+NC=N
extern double V,L; // box volume and size V=L^3;
extern double sigmaA,sigmaB,sigmaC;
extern double sigmaAB,sigmaBC,sigmaCA;
extern double rc,rc2; // rcutoff and rcutoff^2

extern double RX[MAX_N],RY[MAX_N],RZ[MAX_N]; // particle coordinates 
extern double sigma[MAX_N]; // particle diameter
extern int identity[MAX_N]; // particle identity 1:A 2:B 3:C

extern int NNeighbor[MAX_N]; // number of neighbors of each particle
extern int NeighborList[MAX_N][MAX_Neighbors]; // neighbor list of each particle

extern double dRX[MAX_Neighbors],dRY[MAX_Neighbors],dRZ[MAX_Neighbors]; // displacement of neighbor from one particle
extern double dRij[MAX_Neighbors],dRij2[MAX_Neighbors];
extern int Tag[MAX_Neighbors]; // identity of the neighbor of a center particle
extern int Edges[MAX_Neighbors]; // number of edges of the face between a neighbor and a center particle

extern int SoluteCheck;
extern int NF3,NF4,NF5,NF6,NF7,NF8,NF9;
extern int NFaces; // number of faces of a center particle
extern int NEdges; //  number of edges of a center particle
extern int NVertices;//number of vertices of a center particle
extern int IVer[MAX_Vertices],JVer[MAX_Vertices],KVer[MAX_Vertices];// three particles indices that enclose a vertex
extern double dRVX[MAX_Vertices],dRVY[MAX_Vertices],dRVZ[MAX_Vertices];// relative displacement of vertex from the center particle

void ReadInput(void); // read particle number, size, coordinates etc
void Writemovie(FILE *FilePtr); // write .xyz file
void BubbleSort(int ncan);
void Work(int I0,int ncan);
int ShareEdge(int i,int j); // whether two vertices connected edge
