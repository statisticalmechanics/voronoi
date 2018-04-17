#include "system.h"

int Nsample,Nstart,Nend,samplefrequency;

int JobIndex;
int N,NA,NB,NC; // number of particles in ternary mixtures: A,B,C. NA+NB+NC=N
double V,L; // box volume and size V=L^3;
double sigmaA,sigmaB,sigmaC;
double sigmaAB,sigmaBC,sigmaCA;
double rc,rc2; // rcutoff and rcutoff^2

double RX[MAX_N],RY[MAX_N],RZ[MAX_N]; // particle coordinates 
double sigma[MAX_N]; // particle diameter
int identity[MAX_N]; // particle identity 1:A 2:B 3:C

int NNeighbor[MAX_N]; // number of neighbors of each particle
int NeighborList[MAX_N][MAX_Neighbors]; // neighbor list of each particle

double dRX[MAX_Neighbors],dRY[MAX_Neighbors],dRZ[MAX_Neighbors]; // displacement of neighbor from one particle
double dRij[MAX_Neighbors],dRij2[MAX_Neighbors];
int Tag[MAX_Neighbors]; // identity of the neighbor of a center particle
int Edges[MAX_Neighbors]; // number of edges of the face between a neighbor and a center particle

int SoluteCheck;
int NF3,NF4,NF5,NF6,NF7,NF8,NF9;
int NFaces; // number of faces of a center particle
int NEdges; //  number of edges of a center particle
int NVertices;//number of vertices of a center particle
int IVer[MAX_Vertices],JVer[MAX_Vertices],KVer[MAX_Vertices];// three particles indices that enclose a vertex
double dRVX[MAX_Vertices],dRVY[MAX_Vertices],dRVZ[MAX_Vertices];// relative displacement of vertex from the center particle
