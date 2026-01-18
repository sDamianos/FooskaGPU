#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"
#include "bc.h"
#include "cons2prim.h"


__global__ void convert2PrimVecProx(INPUT *d_input, int start, int nFcs, int *fc2el, double *consVec, double *primVecF, int max, int *neigRank4fc, int *lclProx4fc, double *sendBuff){
	
	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;

	double consVecReg[10];
	double sol[10];
	int neigRank, iProx;
	
	while(ifc < nFcs){
	
		for(int iVar = 0; iVar < 10; iVar++){
			consVecReg[iVar] = consVec[fc2el[ifc*2+0]*10 + iVar];
		}			
		cons2prim(d_input,consVecReg,sol);
		
		for(int iVar = 0; iVar < 10; iVar++){
			primVecF[ifc*10*2 + 0*10 + iVar] = sol[iVar];
		}
			
		neigRank = neigRank4fc[ifc];
		iProx    = lclProx4fc[ifc];
		for(int iVar = 0; iVar < 10; iVar++){
			sendBuff[neigRank*max*10 + iProx*10 + iVar] = sol[iVar];
		}

		ifc += blockDim.x * gridDim.x;
	}
}


__global__ void convert2PrimVec(INPUT *d_input, int start, int nFcs, int *fc2el, int *boundCond,  float *n, float * nt1, float *nt2, double *consVec, double *primVecF){
	
	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;
	int id =0, bc = 0;

	double consVecReg[10];
	float nVec[3],nt1Vec[3],nt2Vec[3];
	double sol[10];
	
	while(ifc < nFcs){
		
		bc      = boundCond[ifc];
		id      = (bc == 0);

		for(int iVar = 0; iVar < 10; iVar++){
			consVecReg[iVar] = consVec[fc2el[ifc*2+0]*10 + iVar];
		}	
		cons2prim(d_input,consVecReg,sol);
		//if(fc2el[ifc*2+0]==17022||fc2el[ifc*2+1]==17022){printf("ifc %d iel %d avfD %.10le sol[0] %.10le are1 are2 %f %f\n", ifc, fc2el[ifc*2+0], consVecReg[0], sol[6], sol[7]);}
		
		for(int iVar = 0; iVar < 10; iVar++){
			primVecF[ifc*10*2 + 0*10 + iVar] = sol[iVar];
		}
			
		for(int iVar = 0; iVar < 10; iVar++){
			consVecReg[iVar] = consVec[fc2el[ifc*2+id]*10 + iVar];
		}
		cons2prim(d_input,consVecReg,sol);
		//if(fc2el[ifc*2+0]==17022||fc2el[ifc*2+1]==17022){printf("ifc %d iel %d  sol[0] %.10le are1 are2 %f %f\n", ifc,fc2el[ifc*2+id], consVecReg[0], sol[6], sol[7]);}
	
		for(int idr =0; idr < 3; idr++){
			nVec[idr] = n[ifc*3 + idr];
			nt1Vec[idr] = nt1[ifc*3 + idr];
			nt2Vec[idr] = nt2[ifc*3 + idr];
		}

		if(bc != 0){bcValues(d_input,nVec,nt1Vec,nt2Vec,sol,bc);}
		
		for(int iVar = 0; iVar < 10; iVar++){
			primVecF[ifc*10*2 + 1*10 + iVar] = sol[iVar];
		}
		
		ifc += blockDim.x * gridDim.x;
	}
}

void computePrimVec(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);
		
	convert2PrimVecProx<<<(nProxTh+511)/512,512,0,sProx>>>(d_input, nInner, nFcs, mesh->fc2el,field->consVec,field->primVecF, comm->nProxFacesMax, comm->neigRank4fc,comm->lclProx4fc, comm->sendbuff);
	
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);	
//	cudaEventRecord(start,sMain);	

	convert2PrimVec<<<(nInner+511)/512,512,0,sMain>>>(d_input, 0 , nInner, mesh->fc2el,mesh->boundCond,mesh->n,mesh->nt1,mesh->nt2,field->consVec,field->primVecF); 

//	cudaEventRecord(stop, sMain);
//	cudaEventSynchronize(stop);
//	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("PrimVec time %3.5f ms\n", elapsedTime);
	
	cudaStreamSynchronize(sProx);
	sendArray(comm,0,input.NEQ,comm->sendbuff);		
}

