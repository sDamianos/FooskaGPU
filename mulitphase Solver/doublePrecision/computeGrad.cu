#include "strdata.h"
	
__global__ void initializeLimiters(int nElem, double *wMin, double *wMax){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while(i < nElem*10){
		wMin[i] = 1e12;
		wMax[i] = -1e12;

		i += blockDim.x * gridDim.x;
	}
}
	


__global__ void initializeGrad(int nElem, double *grad){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while( i < nElem*10*3){
		grad[i] = 0.0;

		i += blockDim.x * gridDim.x;
	}

}

__device__ void atomicMinMaxDouble(double* addrMin, double* addrMax, double value) {
	// Update maximum
	unsigned long long int* addrMax_ull = (unsigned long long int*)addrMax;
	unsigned long long int oldMax = *addrMax_ull, assumedMax;
	do {
	assumedMax = oldMax;
		double oldValMax = __longlong_as_double(assumedMax);
		if (oldValMax >= value) break; // Already larger or equal
		oldMax = atomicCAS(addrMax_ull, assumedMax, __double_as_longlong(value));
	} while (assumedMax != oldMax);

	// Update minimum
	unsigned long long int* addrMin_ull = (unsigned long long int*)addrMin;
	unsigned long long int oldMin = *addrMin_ull, assumedMin;
	do {
		assumedMin = oldMin;
		double oldValMin = __longlong_as_double(assumedMin);
		if (oldValMin <= value) break; // Already smaller or equal
		oldMin = atomicCAS(addrMin_ull, assumedMin, __double_as_longlong(value));
	} while (assumedMin != oldMin);
}

__global__ void greenGauss(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, double *n, double *area, double *volume, double *primVecF, double *grad, double *wMin, double *wMax){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/10;
	int iv  = tid%10;
	int bc;
	
	double areaF;
	double vol[2];

	double nVec[3];
	double solL, solR;

	while(ifc < nFcs){
		
		solL = primVecF[ifc*10*2 + 0*10 + iv]; 
		solR = primVecF[ifc*10*2 + 1*10 + iv]; 
		
		nVec[0] = n[ifc*3 + 0]; 
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];
		
		bc = boundCond[ifc];

		atomicMinMaxDouble(&wMin[fc2el[ifc*2+0]*10 + iv], &wMax[fc2el[ifc*2+0]*10 + iv], solR);
		if(bc==0){atomicMinMaxDouble(&wMin[fc2el[ifc*2+1]*10 + iv], &wMax[fc2el[ifc*2+1]*10 + iv], solL);}

		areaF = area[ifc];	
		vol[0] = volume[fc2el[ifc*2+0]];
		if(bc==0){vol[1] = volume[fc2el[ifc*2+1]];}
		
		for(int idr = 0; idr < 3; idr++){
			atomicAdd( &grad[fc2el[ifc*2+0]*10*3+iv*3 + idr], ((solL+solR)/2)*nVec[idr]*areaF/vol[0]);
			if(bc==0){atomicAdd( &grad[fc2el[ifc*2+1]*10*3+iv*3 + idr], - ((solL+solR)/2)*nVec[idr]*areaF/vol[1]);}
		}
		
		tid += blockDim.x * gridDim.x;
		ifc = tid/10;
		iv  = tid%10;
	}
}


__global__ void greenGaussProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, double *n, double *area, double *volume, int *neigRank4fc, int *lclFc2idRcv, double *primVecF, double *grad, double *wMin, double *wMax, double *recvBuff){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = start + tid/10;
	int iv  = tid%10;
	
	double areaF;
	double vol[2];

	double nVec[3];
	double solL, solR;
	int neigRank, iProx;

	while(ifc < nFcs){
	
		solL = primVecF[ifc*10*2 + 0*10 + iv]; 
		neigRank = neigRank4fc[ifc];
		iProx    = lclFc2idRcv[ifc];
		solR = recvBuff[neigRank*max*10 + iProx*10 + iv];
	
		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];
		
		atomicMinMaxDouble(&wMin[fc2el[ifc*2+0]*10 + iv], &wMax[fc2el[ifc*2+0]*10 + iv], solR);

		areaF = area[ifc];	
		vol[0] = volume[fc2el[ifc*2+0]];
		
		for(int idr = 0; idr < 3; idr++){
			atomicAdd( &grad[fc2el[ifc*2+0]*10*3+iv*3 + idr], ((solL+solR)/2)*nVec[idr]*areaF/vol[0]);
		}
		
		tid += blockDim.x * gridDim.x;
		ifc += tid/10;
		iv  = tid%10;
	}
}

__global__ void collectBuff(int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclProx4Fc, double *grad, double *sendBuff){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;	
	int ifc = start + tid/(3*3);
	int iv  = (tid-ifc)/3; 
	int idr = (tid-ifc)%3;     


	int neigRank, iProx;
	while(ifc < nFcs){
		neigRank = neigRank4fc[ifc];
		iProx    = lclProx4Fc[ifc];
		
		sendBuff[neigRank*max*3*3 + iProx*3*3 + iv*3 + idr] = grad[fc2el[ifc*2+0]*10*3+ (iv+3)*3 + idr];
		
		tid += blockDim.x * gridDim.x;
		ifc += tid/(3*3);
		iv   = (tid-ifc)/3;
		idr  = (tid-ifc)%3;
	}
}


void computeGrad(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitGrads, cudaEvent_t evGreenGauss){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);
	
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);	
//	cudaEventRecord(start,sMain);	

	initializeLimiters<<<(mesh->nElem*10+511)/512,512,0,sMain>>>(mesh->nElem, field->wMin, field->wMax); 

	initializeGrad<<<(mesh->nElem*10*3+511)/512,512,0,sMain>>>(mesh->nElem,field->grad); 
	cudaEventRecord(evInitGrads, sMain);
	
	greenGauss<<<(nInner*10+511)/512,512,0,sMain>>>(d_input, nInner, mesh->fc2el, mesh->boundCond, mesh->n, mesh->area, mesh->volume,field->primVecF, field->grad, field->wMin, field->wMax); 

//	cudaEventRecord(stop, sMain);
//	cudaEventSynchronize(stop);
 //	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("Grad time %3.5f ms\n", elapsedTime);

	MPI_Waitall(comm->nNeigRanks, comm->recvRequests[0], MPI_STATUSES_IGNORE);	
	cudaStreamWaitEvent(sProx, evInitGrads, 0);
	cudaStreamSynchronize(sProx);
	greenGaussProx<<<(nProxTh*10+511)/512,512,0,sProx>>>(d_input, nInner, nFcs, comm->nProxFacesMax, mesh->fc2el, mesh->n, mesh->area, mesh->volume, comm->neigRank4fc, comm->lclFc2idRcv, field->primVecF, field->grad, field->wMin, field->wMax, comm->recvbuff);	
	cudaEventRecord(evGreenGauss, sProx);
			
	collectBuff<<<(nProxTh*3*3+511)/512,512,0,sProx>>>(nInner, nFcs, comm->nProxFacesMax, mesh->fc2el, comm->neigRank4fc, comm->lclProx4fc, field->grad, comm->sendBuffGrad);
	
	cudaStreamSynchronize(sProx);
	sendArray(comm,2,3*3,comm->sendBuffGrad);

}

