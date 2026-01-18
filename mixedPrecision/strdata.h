#include <cuda_runtime.h>     
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
	((double *)(A))[3*0+0], ((double *)(A))[3*1+0], ((double *)(A))[3*2+0],		\
	((double *)(A))[3*0+1], ((double *)(A))[3*1+1], ((double *)(A))[3*2+1], 	\
	((double *)(A))[3*0+2], ((double *)(A))[3*1+2], ((double *)(A))[3*2+2]) )

extern int iStep;
extern int iStart;
extern double physicalTime;


struct INPUT{
	float dt;
	float cfl;
	int adaptDt;
	int adaptDtItr;
	float maxDt;
	double avfMin;
	
	int rkSteps;
	float finalTime;
	float exportTime;
	int restartStep;
	int restartFlag;
	int neglectX;
	int neglectY;
	int neglectZ;
	int nMaterials;
	int NEQ;
	double scale;
	int viscFlag;
	int turbulence;
	int woodsFlag;
	
	int thermFntRlxFlag;
	float thermRlxRate;
	float thermUrf;
	int thermItrMax;
	
	int chemFlag;
	int chemFntRlxFlag;
	float chemRlxRate;	
	float chemUrf;
	int chemItrMax;
	int chemCondesation;
	float chemLimit;

	int consCalcFlag;
	int consCalcIvar;

	float satTempUrf;

	float pInf[2];
	double gamma[2];
	float Cp[2];
	float eta[2];
	float satConst[2];
	float dynVisc[2];
	float heatCoef[2];
	float Pr[2];
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
}MESH;

typedef struct{
	double *consVec0;
	double *consVec;
	double *primVecF;
	float *RHS;
	float *grad;
	double *wMax;
	double *wMin;
	float *extensionVec;
	float *theta;
	int *idTherm2iel;
	int *idChem2iel;
	int *nTherm;
	int *nChem;
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
__host__ void translator_ptw(double scale, int *nElem, int *nGlbNodes, int **nFaces, int **nLclNodes, int **lcl2GlbNd, int **boundCond, double **cellCenter, double **fcCenter, double **glbNdCrdnts);
__host__ void meshConnectivity(int nElem,int nGlbNodes,int *nFaces,int *nLclNodes,int *lcl2GlbNd, int *neig, int *neigFc, int *boundCond, int *nEl4GlbNd, int **idGlbEl4GlbNd, int *maxEl4GlbNd, double *fcCenter);
__host__ void allocateHost(int nElem,int nGlbNodes, int **neig, int **h_neigFc, int **h_nEl4GlbNd, int **h_part);
__host__ void decompose(int nElem, int *nFaces, int *d_neig, int *d_part);

__host__ void buildLocalMesh(HOSTLOCALDATA *hLclData, int nElemGlb, int nGlbNodes, int max, int *h_nLclNodes, int *h_neig, int *h_neigFc, int *h_nFaces, int *h_part, int *h_boundCond, int *h_nEl4GlbNd, int *h_idGlbEl4GlbNd,  int *h_elNd2GlbNd, double *h_cellCenter, double *h_fcCenter, double *h_glbNdCord);

__host__ void writeMeshRestart(HOSTLOCALDATA *hLclData, int rank, int size, int nGlbElem, int nGlbNodes);
__host__ void writeFieldRestart(int rank, MESH *mesh, FIELD *field);
__host__ void restartDomain(MESH *mesh, EXPORT *exp, COMM *comm, FIELD *field);
__host__ void readMesh(HOSTLOCALDATA *hLclData, int rank);
__host__ void readComm(HOSTLOCALDATA *hLclData, int rank, int size);
__host__ void readField(int nElem, int rank, double *d_consVec);
__host__ void communicateVolumes(HOSTLOCALDATA *hLclData, int rank, int size, int nNeigRanks, int max, int *neigRanks, int *nProxFaces, int *neigRank4fc, int *lclProx4fc, int * lclFc2idRcv, int *fc2el, float *volume);

__host__ void partitionConnectivity(int nLclElem, int rank, int size, int nNeigRanks, int max, int *nFaces, int *lcl2glbEl, int *lcl2glbFc, int *glb2lclFc, int *neig, int *part, int *elFc2fc, int *neigRanks, int *nProxFaces, int *neigRank4fc, int *lclProx4fc, int *lclFc2idRcv);

__host__ void meshProperties(HOSTLOCALDATA *hLclData, int *nLclNodes, int *fc2el, int *fc2elFc, int *el2glbNd);
__host__ void startDevice(int nGlbElem, int nGlbNodes, int rank, int size, HOSTLOCALDATA *hLclData, MESH *mesh, FIELD *field, COMM *comm, EXPORT *exp);
__host__ void initialize(int nElem, float *cellCenter, double *consVec);
__host__ void export_plt(EXPORT *exp, FIELD *field);
__host__ void explicitLoop(MESH *mesh, EXPORT *exp, COMM *comm, FIELD *field);

__host__ void computePrimVec(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx);
__host__ void computeGrad(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitGrads, cudaEvent_t evGreenGauss);
__host__ void computeFaceValues(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evGreenGauss, cudaEvent_t evCompLimiter);
__host__ void computeRHS(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitRHS, cudaEvent_t evFluxDone);
__host__ void initRk(MESH *mesh, FIELD *field, cudaStream_t sMain);
__host__ void update(int irk, MESH *mesh, FIELD *field, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evFluxDone, cudaEvent_t evStepDone);
__host__ void relaxation(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evStepDone);
__host__ void computeTimeStep(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx);

__global__ void computeViscousStresses(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, double *primVecF, float *grad, float *n, float *area, float *volume, float *RHS);
__global__ void computeViscousStressesProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclFc2idRcv, double *primVecF, float *grad, float *n, float *area, float *volume, float *RHS, float *recvBuffGrad, float *recvBuff);
__global__ void thermalRlx(INPUT *d_input, int nElem, int itrMax, float urf, int *idTherm2iel, double *consVec);
__global__ void chemRlx(INPUT *d_input, int nElem, int itrMax, float urf, int *idiChem2iel, double *consVec);

__host__ void sendArray(COMM * comm, int iField, int size, float *arr);
__host__ void recvArray(COMM * comm, int iField, int size, float *arr);
__host__ void consCalc(MESH *mesh, FIELD *field, cudaStream_t sMain);
