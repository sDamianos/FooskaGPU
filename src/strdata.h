#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(name, mode) _mkdir(name)
#endif

#include"metis.h"
#include"parmetis.h"

/* Define macros */
#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


#ifndef sign
#define sign(a,b) ((b)<0 ? -fabs(a) : fabs(a))
#endif

#define printmat3_fmt(A)	#A " =\n%g %g %g\n%g %g %g\n%g %g %g\n"

#define printmat3(A)	( printf(printmat3_fmt(A), 					\
	((float *)(A))[3*0+0], ((float *)(A))[3*1+0], ((float *)(A))[3*2+0],		\
	((float *)(A))[3*0+1], ((float *)(A))[3*1+1], ((float *)(A))[3*2+1], 	\
	((float *)(A))[3*0+2], ((float *)(A))[3*1+2], ((float *)(A))[3*2+2]) )

extern int iStep;
extern int iStart;
extern float physicalTime;


struct INPUT{
	int turbFlag;
	int rkSteps;
	int nSteps;
	int exportStep;
	int restartStep;
	int restartFlag;
	int neglectX;
	int neglectY;
	int neglectZ;
	int nMaterials;
	int NEQ;
	int turbulence;
	float scale;
	float dt;
	float  pInf;
	float  gamma;
	float  dynVisc;
	float  heatCoef;
	float  diffCoef;
	float  Cp;
	float  Pr;
};

typedef struct{
	int   nElem;      
	int   nGhostElem;
	int   nFcs;       
	int   nOuter;
	int   nProxTot;
	int   *fc2el;     
	int   *boundCond;  
	float *cellCenter; 
	float *fcCenter;   
	float *n;          
	float *nt1;        
	float *nt2;        
	float *area;       
	float *volume;     
} MESH;

typedef struct{
	float *consVec0;
	float *consVec;
	float *primVecF;
	float *RHS;
	float *grad;
	float *wMax;
	float *wMin;
	float *extensionVec;
	float *theta;
} FIELD;

typedef struct{	
	int nNeigRanks;    	     
	int size;	   	   
	int rank;          	  
	int nProxFacesMax; 	    
	int *neigRank4fc;
	int *lclProx4fc;
	int *lclFc2idRcv;
	int *nProxFaces;	    
	int *neigRanks;		    
	float *sendbuff;	  
	float *sendBuffGrad;
	float *recvbuff;
	float *recvBuffGrad;	  
	MPI_Request **recvRequests;
}COMM;

typedef struct{
	int nGlbNodes;
	int nGlbElem;
	int nNodes;       
	int nElem;
	int maxEl4nd;  
	int *lcl2glbEl;
	int *lcl2glbNd;  
	int *elNd2lclNd;
	int *nElNds;
	int *nEl4nd;    
	int *idEl4nd;     
	float *glbNdCord; 
}EXPORT;

typedef struct{  
	int    nElem;      
	int    nGhostElem;
	int    nOuter;
	int    nProxTot;
	int    nFcs;       
	int   *fc2el;      
	int   *boundCond;  
	float *cellCenter; 
	float *fcCenter;   
	float *n;          
	float *nt1;        
	float *nt2;        
	float *area;       
	float *volume;     
	
	int nNeigRanks;    	     
	int nProxFacesMax; 	    
	int *neigRank4fc;
	int *lclProx4fc;
	int *lclFc2idRcv;
	int *nProxFaces;	    
	int *neigRanks;		    
	MPI_Request **recvRequests;  
	
	int nNodes;       
	int maxEl4nd;  
	int *lcl2glbEl;
	int *lcl2glbNd;  
	int *nElNds;
	int *elNd2lclNd;
	int *nEl4nd;     
	int *idEl4nd;     
	float *glbNdCord; 
}HOSTLOCALDATA;


extern struct INPUT input;
extern struct INPUT *d_input;


__host__ void buildInput();
__host__ void mainStart(MESH *mesh,EXPORT *exp, COMM *comm, FIELD *field);
__host__ void translator_ptw(float scale, int *nElem, int *nGlbNodes, int **nFaces, int **nLclNodes, int **lcl2GlbNd, int **boundCond, float **cellCenter, float **fcCenter, float **glbNdCrdnts);
__host__ void meshConnectivity(int nElem,int nGlbNodes,int *nFaces,int *nLclNodes,int *lcl2GlbNd, int *neig, int *neigFc, int *boundCond, int *nEl4GlbNd, int **idGlbEl4GlbNd, int *maxEl4GlbNd, float *fcCenter);
__host__ void allocateHost(int nElem,int nGlbNodes, int **neig, int **h_neigFc, int **h_nEl4GlbNd, int **h_part);
__host__ void decompose(int nElem, int *nFaces, int *d_neig, int *d_part);

__host__ void buildLocalMesh(HOSTLOCALDATA *hLclData, int nElemGlb, int nGlbNodes, int max, int *h_nLclNodes, int *h_neig, int *h_neigFc, int *h_nFaces, int *h_part, int *h_boundCond, int *h_nEl4GlbNd, int *h_idGlbEl4GlbNd,  int *h_elNd2GlbNd, float *h_cellCenter, float *h_fcCenter, float *h_glbNdCord);

__host__ void writeMeshRestart(HOSTLOCALDATA *hLclData, int rank, int size, int nGlbElem, int nGlbNodes);
__host__ void writeFieldRestart(int rank, MESH *mesh, FIELD *field);
__host__ void restartDomain(MESH *mesh, EXPORT *exp, COMM *comm, FIELD *field);
__host__ void readMesh(HOSTLOCALDATA *hLclData, int rank);
__host__ void readComm(HOSTLOCALDATA *hLclData, int rank, int size);
__host__ void communicateVolumes(HOSTLOCALDATA *hLclData, int rank, int size, int nNeigRanks, int max, int *neigRanks, int *nProxFaces, int *neigRank4fc, int *lclProx4fc, int * lclFc2idRcv, int *fc2el, float *volume); 
__host__ void readField(int nElem, int rank, float *d_consVec);

__host__ void partitionConnectivity(int nLclElem, int rank, int size, int nNeigRanks, int max, int *nFaces, int *lcl2glbEl, int *lcl2glbFc, int *glb2lclFc, int *neig, int *part, int *elFc2fc, int *neigRanks, int *nProxFaces, int *neigRank4fc, int *lclProx4fc, int *lclFc2idRcv);

__host__ void meshProperties(HOSTLOCALDATA *hLclData, int *nLclNodes, int *fc2el, int *fc2elFc, int *el2glbNd);
__host__ void startDevice(int nGlbElem, int nGlbNodes, int rank, int size, HOSTLOCALDATA *hLclData, MESH *mesh, FIELD *field, COMM *comm, EXPORT *exp);
__host__ void initialize(int nElem, float *cellCenter, float *consVec);
__host__ void export_plt(EXPORT *exp, float *consVec);
__host__ void explicitLoop(MESH *mesh, EXPORT *exp, COMM *comm, FIELD *field);

__host__ void computePrimVec(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx);
__host__ void computeGrad(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitGrads, cudaEvent_t evGreenGauss);
__host__ void computeFaceValues(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evGreenGauss, cudaEvent_t evCompLimiter);
__host__ void computeRHS(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitRHS, cudaEvent_t evFluxDone);
__host__ void initRk(MESH *mesh, FIELD *field, cudaStream_t sMain);
__host__ void update(int irk, MESH *mesh, FIELD *field, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evFluxDone, cudaEvent_t evStepDone);



__global__ void computeViscousStresses(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, float *primVecF, float *grad, float *n, float *area, float *volume, float *RHS);
__global__ void computeViscousStressesProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclFc2idRcv, float *primVecF, float *grad, float *n, float *area, float *volume, float *RHS, float *recvBuffGrad, float *recvBuff);

__host__ void sendArray(COMM * comm, int iField, int size, float *arr);
__host__ void recvArray(COMM * comm, int iField, int size, float *arr);
