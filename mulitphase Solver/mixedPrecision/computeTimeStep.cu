#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"

__global__ void initializeCellSignal(int nElem, float *theta){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while(i < nElem){
		theta[i] = -1e9;
		i += blockDim.x * gridDim.x;
	}
}

__device__ float computeFaceMaxSignal(INPUT *d_input, double *primVec, float n[3], float nt1[3], float nt2[3]){

	float rhoLiquid, rhoGas, pLiquid, pGas;
	float cLiquid, cGas, cMax;
	float velGlb[3], velLcl[3];
	float cMix, avf, rhoMix, YL, YG;
	
	avf	  = primVec[0];
	rhoLiquid = primVec[1];
	rhoGas	  = primVec[2];
	velGlb[0] = primVec[3];
	velGlb[1] = primVec[4];
	velGlb[2] = primVec[5];
	pLiquid   = primVec[6];
	pGas      = primVec[7];
	rhoMix    = avf*rhoLiquid + (1-avf)*rhoGas;
	cLiquid  = eqos(d_input,0,4,rhoLiquid, pLiquid);
	cGas     = eqos(d_input,1,4,rhoGas, pGas);
	
	if(primVec[8] >  d_input->chemLimit && d_input->woodsFlag == 1){
		cMix = 1/sqrt(rhoMix*(avf/(rhoLiquid*cLiquid*cLiquid) + (1-avf)/(rhoGas*cGas*cGas)));
	}else{	
		YL = avf*rhoLiquid/rhoMix;
		YG = (1-avf)*rhoGas/rhoMix;
		cMix = sqrt(YL*cLiquid*cLiquid + YG*cGas*cGas);
	}
	
	rotateVector_glb2lcl(velGlb,velLcl,n,nt1,nt2);
//	cMax     = fmax(fabs(velLcl[0])+cLiquid, fabs(velLcl[0]) + cGas);
	cMax     = fabs(velLcl[0])+cMix;
	
//	printf("rhoLiq %f pLiq %f us %f c %f\n", rhoLiquid, pLiquid, velLcl[0], cLiquid);

	return cMax;
}

__device__ float atomicMaxFloat(float* addr, float value) { 
	int* addr_i = (int*)addr; 
	int old = *addr_i, assumed; 
	do { 
		assumed = old; 
		float old_val = __int_as_float(assumed); 
		if (old_val >= value) break; 
		old = atomicCAS(addr_i, assumed, __float_as_int(value)); 
	} while (assumed != old); 
	return __int_as_float(old); 
}

/*
__global__ void computeLocalSignalSumProx(INPUT *d_input, int nInner, int nFcs, int *fc2el, float *n, float *nt1, float *nt2, float *primVecF, float *area, float *cellSignalSum, int max, int *neigRank4fc, int *lclFc2idRcv, float *recvBuff){

	int ifc = nInner + threadIdx.x + blockIdx.x * blockDim.x;
	
	float cMax, cLeftMax, cRightMax, nVec[3], nt1Vec[3], nt2Vec[3];
	int neigRank, iProx;
	float areaF;

	cLeftMax = cRightMax = 0;
	while(ifc < nFcs){
				
		for(int idr = 0; idr < 3; idr++){
			nVec[idr] = n[ifc*3 + idr];
			nt1Vec[idr] = nt1[ifc*3 + idr];
			nt2Vec[idr] = nt2[ifc*3 + idr];
		}	
		
		cLeftMax = computeFaceMaxSignal(d_input,&primVecF[ifc*2*10 + 0*10],nVec,nt1Vec,nt2Vec);

		neigRank = neigRank4fc[ifc];
		iProx    = lclFc2idRcv[ifc];
		cRightMax = computeFaceMaxSignal(d_input,&recvBuff[neigRank*max*10 + iProx*10],nVec,nt1Vec,nt2Vec);
		
		cMax = fmax(cLeftMax, cRightMax);	

		areaF = area[ifc];
		atomicAdd(&cellSignalSum[fc2el[ifc*2+0]], cMax*areaF);

		ifc += blockDim.x*gridDim.x;
	}
	
}
*/

__global__ void computeLocalSignalSum(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, float *n, float *nt1, float *nt2, double *primVecF, float *area, float *cellSignalSum){

	int ifc = threadIdx.x + blockIdx.x * blockDim.x;
	
	float cMax, cLeftMax, cRightMax, nVec[3], nt1Vec[3], nt2Vec[3];
	int bc;
//	float areaF;

	cLeftMax = cRightMax = 0;
	while(ifc < nFcs){
				
		for(int idr = 0; idr < 3; idr++){
			nVec[idr] = n[ifc*3 + idr];
			nt1Vec[idr] = nt1[ifc*3 + idr];
			nt2Vec[idr] = nt2[ifc*3 + idr];
		}	
	
		cLeftMax = computeFaceMaxSignal(d_input,&primVecF[ifc*2*10 + 0*10],nVec,nt1Vec,nt2Vec);

		bc = boundCond[ifc];
		if(bc == 0){
			cRightMax = computeFaceMaxSignal(d_input,&primVecF[ifc*2*10 + 1*10],nVec,nt1Vec,nt2Vec);
		}

		cMax = fmax(cLeftMax, cRightMax);	
//		if(fc2el[ifc*2+0]==0||fc2el[ifc*2+1]==0){printf("c %le\n", cMax);}

		//areaF = area[ifc];
		atomicMaxFloat(&cellSignalSum[fc2el[ifc*2+0]],cMax);
		if(bc==0){atomicMaxFloat(&cellSignalSum[fc2el[ifc*2+1]],cMax);}
//		atomicAdd(           &cellSignalSum[fc2el[ifc*2+0]], cMax*areaF);
//		if(bc==0){atomicAdd( &cellSignalSum[fc2el[ifc*2+1]], cMax*areaF);}

		ifc += blockDim.x*gridDim.x;
	}
	
}

__global__ void computeLocalCFL(int nElem, float *volume, float *cellSignalSum, float *blockMin){

	int tid = threadIdx.x;
	int idx = tid + blockIdx.x * blockDim.x;

	__shared__ float sData[512];
	float localMin = 1e12;

	sData[tid] = 1e9;
	while(idx < nElem){
		float val = cbrt(volume[idx]) / cellSignalSum[idx];
//		if(idx==0){printf("val %le c %le dx %le\n", val, cellSignalSum[idx], pow(volume[idx],0.3333));}
		localMin = fmin(localMin, val);
		idx += blockDim.x * gridDim.x;
	}
	
	sData[tid] = localMin;
	__syncthreads();

	for (int stride = blockDim.x/2; stride > 0; stride >>= 1){
		if (tid < stride){
			sData[tid] = fmin(sData[tid], sData[tid + stride]);
		}
		__syncthreads();
	}

	if (tid == 0){
		blockMin[blockIdx.x] = sData[0];
	}
}


void computeTimeStep(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx){

	int nBlocks = (mesh->nElem+511)/512;
	int nInner  = mesh->nFcs-mesh->nProxTot;
	
	float *cellSignalSum;
	cudaMalloc((void**)&cellSignalSum, mesh->nElem*sizeof(float));

	float *h_blockMin, *d_blockMin;
	h_blockMin = (float*)malloc(nBlocks*sizeof(float));
	cudaMalloc((void**)&d_blockMin, nBlocks*sizeof(float));

	initializeCellSignal<<<(mesh->nElem+511)/512,512,0,sMain>>>(mesh->nElem, cellSignalSum);

	//WAIT STREAM FUNTIONS 
	computeLocalSignalSum<<<(mesh->nFcs+511)/512,512,0,sMain>>>(d_input, nInner, mesh->fc2el, mesh->boundCond, mesh->n, mesh->nt1, mesh->nt2, field->primVecF, mesh->area, cellSignalSum);
	
//	int nProxTh = max(1,mesh->nProxTot);
//	computeLocalSignalSumProx<<<(nProxTh+511)/512,512,0,sMain>>>(d_input, nInner, mesh->nFcs, mesh->fc2el, mesh->n, mesh->nt1, mesh->nt2, field->primVecF, mesh->area, cellSignalSum, comm->nProxFacesMax, comm->neigRank4fc, comm->lclFc2idRcv, comm->recvbuff);

	computeLocalCFL<<<nBlocks,512,0,sMain>>>(mesh->nElem, mesh->volume, cellSignalSum, d_blockMin);
	
	cudaStreamSynchronize(sMain);
	cudaMemcpy(h_blockMin, d_blockMin, nBlocks*sizeof(float), cudaMemcpyDeviceToHost);

	float minGlb;
	float minLcl = 1e15;
	for(int iBlock = 0; iBlock < nBlocks; iBlock++){
		minLcl = fmin(minLcl, h_blockMin[iBlock]);
	//	printf("block %d min %le\n", iBlock, h_blockMin[iBlock]);
	}

	MPI_Allreduce(&minLcl, &minGlb, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

	input.dt = fmin(input.cfl*minGlb, input.maxDt);
	if(comm->rank == 0){printf("Time step set to %le\n", input.dt);}

	cudaFree(cellSignalSum);
	free(h_blockMin);
	cudaFree(d_blockMin);
}
