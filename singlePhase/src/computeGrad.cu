#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"
#include "bc.h"
#include "cons2prim.h"
	
__global__ void initializeLimiters(int nElem, float *wMin, float *wMax){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while(i < nElem*5){
		wMin[i] = 1e12;
		wMax[i] = -1e12;

		i += blockDim.x * gridDim.x;
	}
}
	


__global__ void initializeGrad(int nElem, float *grad){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while( i < nElem*5*3){
		grad[i] = 0.0;

		i += blockDim.x * gridDim.x;
	}

}


__device__ void atomicMinMaxFloat(float* addrMin, float* addrMax, float value) {
	// Update maximum
	int* addrMax_i = (int*)addrMax;
	int oldMax = *addrMax_i, assumedMax;
	do {
		assumedMax = oldMax;
		float oldValMax = __int_as_float(assumedMax);
		if (oldValMax >= value) break;
		oldMax = atomicCAS(addrMax_i, assumedMax, __float_as_int(value));
	} while (assumedMax != oldMax);

	// Update minimum
	int* addrMin_i = (int*)addrMin;
	int oldMin = *addrMin_i, assumedMin;
	do {
		assumedMin = oldMin;
		float oldValMin = __int_as_float(assumedMin);
		if (oldValMin <= value) break;
		oldMin = atomicCAS(addrMin_i, assumedMin, __float_as_int(value));
	} while (assumedMin != oldMin);
}




__global__ void greenGauss(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, float *n, float *area, float *volume, float *primVecF, float *grad, float *wMin, float *wMax){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/5;
	int iv  = tid%5;
	int bc;
	
	float areaF;
	float vol[2];

	float nVec[3];
	float solL, solR;

	while(ifc < nFcs){
		
		solL = primVecF[ifc*5*2 + 0*5 + iv]; 
		solR = primVecF[ifc*5*2 + 1*5 + iv]; 
		
		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];
		
		bc = boundCond[ifc];

		atomicMinMaxFloat(&wMin[fc2el[ifc*2+0]*5 + iv], &wMax[fc2el[ifc*2+0]*5 + iv], solR);
		if(bc==0){atomicMinMaxFloat(&wMin[fc2el[ifc*2+1]*5 + iv], &wMax[fc2el[ifc*2+1]*5 + iv], solL);}

		areaF = area[ifc];	
		vol[0] = volume[fc2el[ifc*2+0]];
		if(bc==0){vol[1] = volume[fc2el[ifc*2+1]];}
		
		for(int idr = 0; idr < 3; idr++){
			atomicAdd( &grad[fc2el[ifc*2+0]*5*3+iv*3 + idr], ((solL+solR)/2)*nVec[idr]*areaF/vol[0]);
			if(bc==0){atomicAdd( &grad[fc2el[ifc*2+1]*5*3+iv*3 + idr], - ((solL+solR)/2)*nVec[idr]*areaF/vol[1]);}
		}
		
		tid += blockDim.x * gridDim.x;
		ifc = tid/5;
		iv  = tid%5;
	}
}

__global__ void greenGaussProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, float *n, float *area, float *volume, int *neigRank4fc, int *lclFc2idRcv, float *primVecF, float *grad, float *wMin, float *wMax, float *recvBuff){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = start + tid/5;
	int iv  = tid%5;
	
	float areaF;
	float vol[2];

	float nVec[3];
	float solL, solR;
	int neigRank, iProx;

	while(ifc < nFcs){
	
		solL = primVecF[ifc*5*2 + 0*5 + iv]; 
		neigRank = neigRank4fc[ifc];
		iProx    = lclFc2idRcv[ifc];
		solR = recvBuff[neigRank*max*5 + iProx*5 + iv];
	
		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];
		
		atomicMinMaxFloat(&wMin[fc2el[ifc*2+0]*5 + iv], &wMax[fc2el[ifc*2+0]*5 + iv], solR);

		areaF = area[ifc];	
		vol[0] = volume[fc2el[ifc*2+0]];
		
		for(int idr = 0; idr < 3; idr++){
			atomicAdd( &grad[fc2el[ifc*2+0]*5*3+iv*3 + idr], ((solL+solR)/2)*nVec[idr]*areaF/vol[0]);
		}
		
		tid += blockDim.x * gridDim.x;
		ifc += tid/5;
		iv  = tid%5;
	}
}

__global__ void collectBuff(int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclProx4Fc, float *grad, float *sendBuff){

	int tid = start*5*3 + threadIdx.x + blockIdx.x * blockDim.x;	
	int ifc = tid/(5*3);
	int iv  = (tid-ifc*5*3)/3; 
	int idr = (tid-ifc*5*3)%3;     


	int neigRank, iProx;
	while(ifc < nFcs){
		neigRank = neigRank4fc[ifc];
		iProx    = lclProx4Fc[ifc];
		
		sendBuff[neigRank*max*5*3 + iProx*5*3 + iv*3 + idr] = grad[fc2el[ifc*2+0]*5*3+iv*3 + idr];
		
		tid += blockDim.x * gridDim.x;
		ifc = tid/(5*3);
		iv  = (tid-ifc*5*3)/3;
		idr = (tid-ifc*5*3)%3;
	}
}

void computeGrad(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitGrads, cudaEvent_t evGreenGauss){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);

	initializeLimiters<<<(mesh->nElem*5+511)/512,512,0,sMain>>>(mesh->nElem, field->wMin, field->wMax); 

	initializeGrad<<<(mesh->nElem*5*3+511)/512,512,0,sMain>>>(mesh->nElem,field->grad); 
	cudaEventRecord(evInitGrads, sMain);
	
	greenGauss<<<(nInner*5+511)/512,512,0,sMain>>>(d_input, nInner, mesh->fc2el, mesh->boundCond, mesh->n, mesh->area, mesh->volume,field->primVecF, field->grad, field->wMin, field->wMax); 
	
	MPI_Waitall(comm->nNeigRanks, comm->recvRequests[0], MPI_STATUSES_IGNORE);	
	cudaStreamWaitEvent(sProx, evInitGrads, 0);
	cudaStreamSynchronize(sProx);
	greenGaussProx<<<(nProxTh*5+511)/512,512,0,sProx>>>(d_input, nInner, nFcs, comm->nProxFacesMax, mesh->fc2el, mesh->n, mesh->area, mesh->volume, comm->neigRank4fc, comm->lclFc2idRcv, field->primVecF, field->grad, field->wMin, field->wMax, comm->recvbuff);	
	cudaEventRecord(evGreenGauss, sProx);
			
	collectBuff<<<(nProxTh*5*3+511)/512,512,0,sProx>>>(nInner, nFcs, comm->nProxFacesMax, mesh->fc2el, comm->neigRank4fc, comm->lclProx4fc, field->grad, comm->sendBuffGrad);
	
	cudaStreamSynchronize(sProx);
	sendArray(comm,2,5*3,comm->sendBuffGrad);
}
