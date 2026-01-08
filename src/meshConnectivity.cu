#include "strdata.h"
#include "operations.h"

__global__ void initializeArray(int *c, int N){
	int ignd = threadIdx.x + blockIdx.x * blockDim.x;
	while (ignd < N) {
		c[ignd] = 0;
		ignd += blockDim.x * gridDim.x;
	}
}


__global__ void compute_nEl4GlbNd(int nGlbElem, int nGlbNodes, int *nLclNodes,int *elNd2GlbNd,int *nEl4GlbNd){
	
	int iel = threadIdx.x + blockIdx.x * blockDim.x;
	int ignd;
	while (iel < nGlbElem){
		for(int ind = 0; ind < nLclNodes[iel]; ind++){
			ignd = elNd2GlbNd[iel*8 + ind];
			atomicAdd( &nEl4GlbNd[ignd], 1 );
		}
		iel += blockDim.x * gridDim.x;
	}
}
	
__global__ void compute_maxEl4glbNd(int nGlbNodes, int *nEl4GlbNd,int *maxPerBlock){

	__shared__ int temp[512];
	
	int ignd = threadIdx.x + blockIdx.x * blockDim.x;
	int i = blockDim.x/2;
	int threadIndex = threadIdx.x;

	temp[threadIndex] = (ignd < nGlbNodes) ? nEl4GlbNd[ignd] : 0;
	__syncthreads();
	
	while(i!=0){
		if(threadIndex < i){
			temp[threadIndex] = max(temp[threadIndex],temp[threadIndex + i]);
		}
		i /= 2;
		__syncthreads();
	}
	
	if(threadIndex == 0){
		maxPerBlock[blockIdx.x] = temp[0];
	}
}
	
__global__ void compute_idGlbEl4GlbNd(int nGlbElem, int maxEl4glbNd, int *elNd2GlbNd,int *nLclNodes,int *nEl4GlbNd,int *idGlbEl4GlbNd){
	
	int iel = threadIdx.x + blockIdx.x * blockDim.x;
	int ignd,slot;
	while (iel < nGlbElem){
		for(int ind = 0; ind < nLclNodes[iel]; ind++){
			ignd = elNd2GlbNd[iel*8 + ind];
			slot = atomicAdd(&nEl4GlbNd[ignd], 1);
			idGlbEl4GlbNd[ignd*maxEl4glbNd + slot] = iel;
		}
		iel += blockDim.x * gridDim.x;
	}
}
	
__global__ void initializeNeigs(int nGlbElem, int *neig_lcl, int *neigFc_lcl){

	int tid =  threadIdx.x + blockIdx.x * blockDim.x;

	while(tid < nGlbElem*6){
		neig_lcl[tid] = 0;
		neigFc_lcl[tid] = 0;

		tid += blockDim.x * gridDim.x;
	}

}

__global__ void	initializeLocalNeig(int nGlbElem, int *neig_lcl, int *neigFc_lcl){

	int tid = threadIdx.x + blockIdx.x *blockDim.x;

	while(tid < nGlbElem*6){
		
		neig_lcl[tid]   = 0;
		neigFc_lcl[tid] = 0;
		
		tid += blockDim.x*gridDim.x;
	}
}
	


// This kernel decrease computational cost by a factor of 13 with the respective part in CPU code. However
// keep in mind that this part is also parallelized in teh CPU version. Due to the face indexing approach
// is used we have contigous memory access, although the non-SoA architecture. If we change to SoA the
// fcCenter arrays although we could have contigous acces when mooving from x coord to y coord we will
// have a huge mmenory jump of nGlbElem*6 which is expect to be much worse than the recent methodology
__global__ void solveConnectivity(int rank, int nGlbElem, int nElemLcl, int maxEl4glbNd, int * nFaces, int * boundCond, int * elNd2GlbNd, int * nLclNodes, int * nEl4GlbNd, int * idGlbEl4GlbNd, float * fcCenter, int * neig, int * neigFc){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = nElemLcl*rank + tid/6;
	int ifc = tid%6; 
	int iel1,ignd;

	while( iel < min(nElemLcl*(rank+1),nGlbElem)){
		if(ifc < nFaces[iel]){
			neig[iel*6 + ifc] = -1;
			if(boundCond[iel*6+ifc] == 0){
				ignd = elNd2GlbNd[iel*8 + d_fcNd2ElNd(ifc,0,nLclNodes[iel])];
				for(int iv = 0; iv < nEl4GlbNd[ignd]; iv++){
					iel1 = idGlbEl4GlbNd[ignd*maxEl4glbNd + iv];
					if(iel1 != iel){
						for(int ifc1 = 0; ifc1 < nFaces[iel1]; ifc1++){
							if( fabs(fcCenter[iel*6*3 + ifc*3 + 0] - fcCenter[iel1*6*3 + ifc1*3 + 0]) < 1e-7 && fabs(fcCenter[iel*6*3 + ifc*3 + 1] - fcCenter[iel1*6*3 + ifc1*3 + 1]) < 1e-7 && fabs(fcCenter[iel*6*3 + ifc*3 + 2] - fcCenter[iel1*6*3 + ifc1*3 + 2]) < 1e-7){
								neig[iel*6+ifc] = iel1;
								neigFc[iel*6+ifc] = ifc1;
								iv   = 100;
								ifc1 = 10;
							}
						}
					}
				}
			}
		}
		tid += blockDim.x * gridDim.x;
		iel  = nElemLcl*rank + tid/6;
		ifc  = tid%6; 
	}
}

void copyToDevice(int nGlbElem,	      int   nGlbNodes,
             	  int   *h_nFaces,    int   **d_nFaces,
                  int   *h_nLclNodes, int   **d_nLclNodes,
                  int   *h_elNd2GlbNd, int   **d_elNd2GlbNd,
		  int   *h_boundCond, int   **d_boundCond,
		  float *h_fcCenter,  float **d_fcCenter,
		  int   **d_neigLcl,  int   **d_neigFcLcl,
		  int   **d_neig,     int   **d_neigFc,
		  int   **d_nEl4GlbNd){
	
	cudaMalloc((void**)d_nFaces,    nGlbElem *         sizeof(int));
	cudaMalloc((void**)d_nLclNodes, nGlbElem *         sizeof(int));
	cudaMalloc((void**)d_elNd2GlbNd, nGlbElem * 8 *     sizeof(int));
	cudaMalloc((void**)d_boundCond, nGlbElem * 6 *     sizeof(int));
	cudaMalloc((void**)d_neigLcl,   nGlbElem * 6 *     sizeof(int)); 
	cudaMalloc((void**)d_neigFcLcl, nGlbElem * 6 *     sizeof(int)); 
	cudaMalloc((void**)d_neig,        nGlbElem * 6 *     sizeof(int)); 
	cudaMalloc((void**)d_neigFc,      nGlbElem * 6 *     sizeof(int)); 
	cudaMalloc((void**)d_nEl4GlbNd,   nGlbNodes *     sizeof(int)); 
	cudaMalloc((void**)d_fcCenter,    nGlbElem * 6 * 3 * sizeof(float));

	
	cudaMemcpy(*d_nFaces,      h_nFaces,      nGlbElem *         sizeof(int),   cudaMemcpyHostToDevice);
	cudaMemcpy(*d_nLclNodes,   h_nLclNodes,   nGlbElem *         sizeof(int),   cudaMemcpyHostToDevice);
	cudaMemcpy(*d_elNd2GlbNd,   h_elNd2GlbNd,   nGlbElem * 8 *     sizeof(int),   cudaMemcpyHostToDevice);
	cudaMemcpy(*d_boundCond,   h_boundCond,   nGlbElem * 6 *     sizeof(int),   cudaMemcpyHostToDevice);
	cudaMemcpy(*d_fcCenter,    h_fcCenter,    nGlbElem * 6 * 3 * sizeof(float), cudaMemcpyHostToDevice);

}
void copyToHost(int nGlbElem,            int nGlbNodes,
		int *h_neig,          int *d_neig,
		int *h_neigFc, 	      int *d_neigFc,
		int *h_nEl4GlbNd,     int *d_nEl4GlbNd,
		int *h_idGlbEl4GlbNd, int *d_idGlbEl4GlbNd,
		int max){
	
	cudaMemcpy(h_neig,          d_neig,          nGlbElem * 6 *       sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_neigFc,        d_neigFc,        nGlbElem * 6 *       sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_nEl4GlbNd,     d_nEl4GlbNd,     nGlbNodes *       sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_idGlbEl4GlbNd, d_idGlbEl4GlbNd, nGlbNodes * max * sizeof(int), cudaMemcpyDeviceToHost);
}

// The GPU architecture for this function decrease by a factor of 10 (!) the computaional cost of this function 
// for the case of 2 Million elements. However the cost of allocating and copying the variables before this functions
// has a non-neglegible cost. As a result the total cost decrease approximatelly by 2, compraing with the CPU architecture
void meshConnectivity(int nGlbElem,int nGlbNodes,int *h_nFaces,int *h_nLclNodes,int *h_elNd2GlbNd, int *h_neig, int *h_neigFc, int *h_boundCond, int *h_nEl4GlbNd, int **h_idGlbEl4GlbNd, int *max, float *h_fcCenter){

	int maxEl4glbNd;
	int rank, size;
	int *d_neigLcl;
	int *d_neigFcLcl;
	int *d_nFaces, *d_nLclNodes, *d_elNd2GlbNd;
	int *d_neig, *d_neigFc, *d_boundCond, *d_nEl4GlbNd;
	int *d_idGlbEl4GlbNd;
	float *d_fcCenter; 
	 	
	copyToDevice(nGlbElem,nGlbNodes,h_nFaces,&d_nFaces,h_nLclNodes,&d_nLclNodes,h_elNd2GlbNd,&d_elNd2GlbNd,h_boundCond,&d_boundCond,h_fcCenter,&d_fcCenter,&d_neigLcl,&d_neigFcLcl,&d_neig,&d_neigFc,&d_nEl4GlbNd);
		
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	initializeArray<<<(nGlbNodes+511)/512,512>>>(d_nEl4GlbNd, nGlbNodes);

	compute_nEl4GlbNd<<<(nGlbElem+511)/512,512>>>(nGlbElem,nGlbNodes,d_nLclNodes,d_elNd2GlbNd,d_nEl4GlbNd);
	
	int * d_maxPerBlock;
	int * h_maxPerBlock;
	cudaMalloc((void**)&d_maxPerBlock, (nGlbNodes+511)/512 * sizeof(int));
	h_maxPerBlock = (int*) malloc((nGlbNodes+511)/512 * sizeof(int));

	compute_maxEl4glbNd<<<(nGlbNodes+511)/512,512>>>(nGlbNodes,d_nEl4GlbNd,d_maxPerBlock);
	cudaMemcpy(h_maxPerBlock,d_maxPerBlock, (nGlbNodes+511)/512 * sizeof(int), cudaMemcpyDeviceToHost);

	maxEl4glbNd = 0;
	for(int iBlk = 0; iBlk < (nGlbNodes+511)/512; iBlk++){
		if(h_maxPerBlock[iBlk] > maxEl4glbNd){
			maxEl4glbNd = h_maxPerBlock[iBlk];
		}
	}
	(*max) = maxEl4glbNd;
	cudaFree(d_maxPerBlock);
	free(h_maxPerBlock);
	
	cudaMalloc((void**)&d_idGlbEl4GlbNd, nGlbNodes*maxEl4glbNd*sizeof(int));
	
	initializeArray<<<(nGlbNodes+511)/512,512>>>(d_nEl4GlbNd, nGlbNodes);

	compute_idGlbEl4GlbNd<<<(nGlbNodes+511)/512,512>>>(nGlbElem,maxEl4glbNd,d_elNd2GlbNd,d_nLclNodes,d_nEl4GlbNd,d_idGlbEl4GlbNd);

//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);	
//	cudaEventRecord(start,0);	

	int nElemLcl = nGlbElem / size + 1;

	initializeLocalNeig<<<(nGlbElem*6+511)/512,512>>>(nGlbElem, d_neigLcl, d_neigFcLcl);
 	
	solveConnectivity<<<(nElemLcl*6+511)/512,512>>>(rank,nGlbElem,nElemLcl,maxEl4glbNd,d_nFaces,d_boundCond,d_elNd2GlbNd,d_nLclNodes,d_nEl4GlbNd,d_idGlbEl4GlbNd,d_fcCenter,d_neigLcl,d_neigFcLcl);
	
	cudaDeviceSynchronize();
	MPI_Barrier(MPI_COMM_WORLD);
		
	MPI_Allreduce(d_neigLcl, d_neig, nGlbElem*6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(d_neigFcLcl, d_neigFc, nGlbElem*6, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		
	cudaFree(d_neigLcl);
	cudaFree(d_neigFcLcl);
	
	*h_idGlbEl4GlbNd = (int*)malloc(nGlbNodes*maxEl4glbNd*sizeof(int));
	copyToHost(nGlbElem, nGlbNodes, h_neig, d_neig, h_neigFc, d_neigFc, h_nEl4GlbNd, d_nEl4GlbNd, *h_idGlbEl4GlbNd, d_idGlbEl4GlbNd, maxEl4glbNd);

	cudaFree(d_nFaces);
	cudaFree(d_nLclNodes);
	cudaFree(d_elNd2GlbNd);
	cudaFree(d_neig);
	cudaFree(d_neigFc);
	cudaFree(d_boundCond);
	cudaFree(d_nEl4GlbNd);
	cudaFree(d_idGlbEl4GlbNd);
	cudaFree(d_fcCenter);

//	cudaEventRecord(stop, 0);
//	cudaEventSynchronize(stop);
//	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("Connectivity Time %3.5f ms\n", elapsedTime);	

}
