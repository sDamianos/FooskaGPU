#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"

__global__ void initializeCellSignal(int nElem, double *theta){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while(i < nElem){
		theta[i] = -1e9;
		i += blockDim.x * gridDim.x;
	}
}

__device__ double computeFaceMaxSignal(INPUT *d_input, double *primVec, double n[3], double nt1[3], double nt2[3]){

	double rhoLiquid, rhoGas, pLiquid, pGas;
	double cLiquid, cGas, cMax;
	double velGlb[3], velLcl[3];
	double cMix, avf, rhoMix, YL, YG;
	
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
	cMax     = fabs(velLcl[0])+cMix;
	
	return cMax;
}

__device__ double atomicMaxDouble(double* addr, double value) {

	unsigned long long int* addr_ull = (unsigned long long int*)addr;
	unsigned long long int old = *addr_ull, assumed;

	do {
		assumed = old;
		double old_val = __longlong_as_double(assumed);
		if (old_val >= value) break;  // Already smaller, no update needed
		old = atomicCAS(addr_ull, assumed, __double_as_longlong(value));
	} while (assumed != old);
	
	return __longlong_as_double(old);  // returns previous value
}


__global__ void computeLocalSignalSumProx(INPUT *d_input, int nInner, int nFcs, int *fc2el, double *n, double *nt1, double *nt2, double *primVecF, double *area, double *cellSignalSum, int max, int *neigRank4fc, int *lclFc2idRcv, double *recvBuff){

	int ifc = nInner + threadIdx.x + blockIdx.x * blockDim.x;
	
	double cMax, cLeftMax, cRightMax, nVec[3], nt1Vec[3], nt2Vec[3];
	int neigRank, iProx;
	double areaF;

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

		atomicMaxDouble(&cellSignalSum[fc2el[ifc*2+0]],cMax);

		ifc += blockDim.x*gridDim.x;
	}
	
}


__global__ void computeLocalSignalSum(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, double *n, double *nt1, double *nt2, double *primVecF, double *area, double *cellSignalSum){

	int ifc = threadIdx.x + blockIdx.x * blockDim.x;
	
	double cMax, cLeftMax, cRightMax, nVec[3], nt1Vec[3], nt2Vec[3];
	int bc;
	double areaF;

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

		atomicMaxDouble(&cellSignalSum[fc2el[ifc*2+0]],cMax);
		if(bc==0){atomicMaxDouble(&cellSignalSum[fc2el[ifc*2+1]],cMax);}

		ifc += blockDim.x*gridDim.x;
	}
	
}

__global__ void computeLocalCFL(int nElem, double *volume, double *cellSignalSum, double *blockMin){

	int tid = threadIdx.x;
	int idx = tid + blockIdx.x * blockDim.x;

	__shared__ double sData[512];
	double localMin = 1e12;

	sData[tid] = 1e9;
	while(idx < nElem){
		double val = pow(volume[idx],0.3333) / cellSignalSum[idx];
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


// COMMENTS: 
// There are two types of CFL computation, 
// (1) the one with the sume of the face weitghted characteristic speed and 
// (2) the other one with the max characteristic sped of the cell and the v^1/3
// here we use the second. THe first one is generaly more cosnervative and is secure
// for high aspect ratio and skewnes off the unstructure mesh. However in good meshes
// reduces significantly the time step much more than is needed (almost one roder of magnitude)

// The speed of sound is not computed thorug the eigen values of the system but
// form the mixture speed of sound used in HLLC.
void computeTimeStep(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx){

	int nBlocks = (mesh->nElem+511)/512;
	int nInner  = mesh->nFcs-mesh->nProxTot;
	
	double *cellSignalSum;
	cudaMalloc((void**)&cellSignalSum, mesh->nElem*sizeof(double));

	double *h_blockMin, *d_blockMin;
	h_blockMin = (double*)malloc(nBlocks*sizeof(double));
	cudaMalloc((void**)&d_blockMin, nBlocks*sizeof(double));

	cudaStreamSynchronize(sMain); //to be sure that data have been received

	initializeCellSignal<<<(mesh->nElem+511)/512,512,0,sMain>>>(mesh->nElem, cellSignalSum);

	computeLocalSignalSum<<<(mesh->nFcs+511)/512,512,0,sMain>>>(d_input, nInner, mesh->fc2el, mesh->boundCond, mesh->n, mesh->nt1, mesh->nt2, field->primVecF, mesh->area, cellSignalSum);
	
	int nProxTh = max(1,mesh->nProxTot);
	computeLocalSignalSumProx<<<(nProxTh+511)/512,512,0,sMain>>>(d_input, nInner, mesh->nFcs, mesh->fc2el, mesh->n, mesh->nt1, mesh->nt2, field->primVecF, mesh->area, cellSignalSum, comm->nProxFacesMax, comm->neigRank4fc, comm->lclFc2idRcv, comm->recvbuff);

	computeLocalCFL<<<nBlocks,512,0,sMain>>>(mesh->nElem, mesh->volume, cellSignalSum, d_blockMin);
	
	cudaStreamSynchronize(sMain);
	cudaMemcpy(h_blockMin, d_blockMin, nBlocks*sizeof(double), cudaMemcpyDeviceToHost);

	double minGlb;
	double minLcl = 1e15;
	for(int iBlock = 0; iBlock < nBlocks; iBlock++){
		minLcl = fmin(minLcl, h_blockMin[iBlock]);
	}

	MPI_Allreduce(&minLcl, &minGlb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	input.dt = fmin(input.cfl*minGlb, input.maxDt);
	if(comm->rank == 0){printf("Time step set to %le\n", input.dt);}

	free(h_blockMin);
	cudaFree(cellSignalSum);
	cudaFree(d_blockMin);
}
