#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"
#include "bc.h"
#include "cons2prim.h"


__global__ void convert2PrimVecProx(INPUT *d_input, int start, int nFcs, int *fc2el, float *consVec, float *primVecF, int max, int *neigRank4fc, int *lclProx4fc, float *sendBuff){
	
	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;

	float consVecReg[5];
	float sol[5];
	int neigRank, iProx;
	
	while(ifc < nFcs){
	
		consVecReg[0] = consVec[fc2el[ifc*2+0]*5 + 0];
		consVecReg[1] = consVec[fc2el[ifc*2+0]*5 + 1];
		consVecReg[2] = consVec[fc2el[ifc*2+0]*5 + 2];
		consVecReg[3] = consVec[fc2el[ifc*2+0]*5 + 3];
		consVecReg[4] = consVec[fc2el[ifc*2+0]*5 + 4];
			
		cons2prim(d_input,consVecReg,sol);
		
		primVecF[ifc*5*2 + 0*5 + 0] = sol[0];
		primVecF[ifc*5*2 + 0*5 + 1] = sol[1];
		primVecF[ifc*5*2 + 0*5 + 2] = sol[2];
		primVecF[ifc*5*2 + 0*5 + 3] = sol[3];
		primVecF[ifc*5*2 + 0*5 + 4] = sol[4];	
			
		neigRank = neigRank4fc[ifc];
		iProx    = lclProx4fc[ifc];
		sendBuff[neigRank*max*5 + iProx*5 + 0] = sol[0];
		sendBuff[neigRank*max*5 + iProx*5 + 1] = sol[1];
		sendBuff[neigRank*max*5 + iProx*5 + 2] = sol[2];
		sendBuff[neigRank*max*5 + iProx*5 + 3] = sol[3];
		sendBuff[neigRank*max*5 + iProx*5 + 4] = sol[4];

		ifc += blockDim.x * gridDim.x;
	}
}


__global__ void convert2PrimVec(INPUT *d_input, int start, int nFcs, int *fc2el, int *boundCond,  float *n, float * nt1, float *nt2, float *consVec, float *primVecF){
	
	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;
	int id =0, bc = 0;

	float consVecReg[5];
	float nVec[3],nt1Vec[3],nt2Vec[3];
	float sol[5];
	
	while(ifc < nFcs){
		
		bc      = boundCond[ifc];
		id      = (bc == 0);
		
		consVecReg[0] = consVec[fc2el[ifc*2+0]*5 + 0];
		consVecReg[1] = consVec[fc2el[ifc*2+0]*5 + 1];
		consVecReg[2] = consVec[fc2el[ifc*2+0]*5 + 2];
		consVecReg[3] = consVec[fc2el[ifc*2+0]*5 + 3];
		consVecReg[4] = consVec[fc2el[ifc*2+0]*5 + 4];
			
		cons2prim(d_input,consVecReg,sol);
		
		primVecF[ifc*5*2 + 0*5 + 0] = sol[0];
		primVecF[ifc*5*2 + 0*5 + 1] = sol[1];
		primVecF[ifc*5*2 + 0*5 + 2] = sol[2];
		primVecF[ifc*5*2 + 0*5 + 3] = sol[3];
		primVecF[ifc*5*2 + 0*5 + 4] = sol[4];	
			
		consVecReg[0] = consVec[fc2el[ifc*2+id]*5 + 0];
		consVecReg[1] = consVec[fc2el[ifc*2+id]*5 + 1];
		consVecReg[2] = consVec[fc2el[ifc*2+id]*5 + 2];
		consVecReg[3] = consVec[fc2el[ifc*2+id]*5 + 3];
		consVecReg[4] = consVec[fc2el[ifc*2+id]*5 + 4];

		cons2prim(d_input,consVecReg,sol);
	
		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];
		
		nt1Vec[0] = nt1[ifc*3 + 0];
		nt1Vec[1] = nt1[ifc*3 + 1];
		nt1Vec[2] = nt1[ifc*3 + 2];

		nt2Vec[0] = nt2[ifc*3 + 0];
		nt2Vec[1] = nt2[ifc*3 + 1];
		nt2Vec[2] = nt2[ifc*3 + 2];
		if(bc != 0){bcValues(d_input,nVec,nt1Vec,nt2Vec,sol,bc);}
		
		primVecF[ifc*5*2 + 1*5 + 0] = sol[0];
		primVecF[ifc*5*2 + 1*5 + 1] = sol[1];
		primVecF[ifc*5*2 + 1*5 + 2] = sol[2];
		primVecF[ifc*5*2 + 1*5 + 3] = sol[3];
		primVecF[ifc*5*2 + 1*5 + 4] = sol[4];

		ifc += blockDim.x * gridDim.x;
	}
}

void computePrimVec(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);
		
	convert2PrimVecProx<<<(nProxTh+511)/512,512,0,sProx>>>(d_input, nInner, nFcs, mesh->fc2el,field->consVec,field->primVecF, comm->nProxFacesMax, comm->neigRank4fc,comm->lclProx4fc, comm->sendbuff);

	convert2PrimVec<<<(nInner+511)/512,512,0,sMain>>>(d_input, 0 , nInner, mesh->fc2el,mesh->boundCond,mesh->n,mesh->nt1,mesh->nt2,field->consVec,field->primVecF); 
	
	cudaStreamSynchronize(sProx);
	sendArray(comm,0,5,comm->sendbuff);
		
}

//	cudaEventRecord(evProxPrimVec, sProx);
//	cudaStreamWaitEvent(0, evProxPrimVec, 0);

