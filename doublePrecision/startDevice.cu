#include "strdata.h"

__global__ void init_ids(int *a, int n) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n){a[i] = i;}
}

__global__ void set_int(int *arr, int val) {
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	if(i == 0){arr[i] = val;}
}


void startDevice(int nGlbElem, int nGlbNodes, int rank, int size, HOSTLOCALDATA *hLclData, MESH *mesh, FIELD *field, COMM *comm, EXPORT *exp){

	int nFcs          = hLclData->nFcs;
	int nElem         = hLclData->nElem;
	int nGhostElem    = hLclData->nGhostElem;
	int nNeigRanks    = hLclData->nNeigRanks;
	int nProxFacesMax = hLclData->nProxFacesMax;
	int nNodes        = hLclData->nNodes;
	int maxEl4nd      = hLclData->maxEl4nd;
	
	// Mesh Struct
	cudaMalloc((void**)&mesh->fc2el,      nFcs * 2 * sizeof(int));
	cudaMalloc((void**)&mesh->boundCond,  nFcs * sizeof(int));
	cudaMalloc((void**)&mesh->cellCenter, nElem * 3 * sizeof(double));
	cudaMalloc((void**)&mesh->fcCenter,   nFcs * 3 * sizeof(double));
	cudaMalloc((void**)&mesh->n,          nFcs * 3 * sizeof(double));
	cudaMalloc((void**)&mesh->nt1,        nFcs * 3 * sizeof(double));
	cudaMalloc((void**)&mesh->nt2,        nFcs * 3 * sizeof(double));
	cudaMalloc((void**)&mesh->area,       nFcs * sizeof(double));
	cudaMalloc((void**)&mesh->volume,     (nElem + nGhostElem) * sizeof(double));
	
	//Field Struct
	cudaMalloc((void**)&field->consVec0,      nElem * input.NEQ * sizeof(double));
	cudaMalloc((void**)&field->consVec,       (nElem + nGhostElem) * input.NEQ * sizeof(double));
	cudaMalloc((void**)&field->primVecF,      nFcs * input.NEQ * 2 * sizeof(double));
	cudaMalloc((void**)&field->RHS,           (nElem + nGhostElem) * input.NEQ * sizeof(double));
	cudaMalloc((void**)&field->grad,          (nElem + nGhostElem) * input.NEQ * 3 * sizeof(double));
	cudaMalloc((void**)&field->wMax,          (nElem + nGhostElem) * input.NEQ * sizeof(double));
	cudaMalloc((void**)&field->wMin,          (nElem + nGhostElem) * input.NEQ * sizeof(double));
	cudaMalloc((void**)&field->extensionVec,  nFcs * input.NEQ * 2 * sizeof(double));
	cudaMalloc((void**)&field->theta,         (nElem + nGhostElem) * input.NEQ * sizeof(double));
	cudaMalloc((void**)&field->idTherm2iel, nElem*sizeof(int));
	init_ids<<<(nElem + 511)/512,512>>>(field->idTherm2iel, nElem);
	cudaMalloc((void**)&field->idChem2iel,  nElem*sizeof(int));
	cudaMalloc(&field->nTherm, sizeof(int));
	set_int<<<1,1>>>(field->nTherm, nElem);
	cudaMalloc(&field->nChem, sizeof(int));

	// Communication Struct
	cudaMalloc((void**)&comm->neigRank4fc,    nFcs * sizeof(int));
	cudaMalloc((void**)&comm->lclProx4fc,     nFcs * sizeof(int));
	cudaMalloc((void**)&comm->lclFc2idRcv,    nFcs * sizeof(int));
	cudaMalloc((void**)&comm->sendbuff,       size * nProxFacesMax * input.NEQ * sizeof(double));
	cudaMalloc((void**)&comm->recvbuff,       size * nProxFacesMax * input.NEQ * sizeof(double));
	cudaMalloc((void**)&comm->sendBuffGrad,   size * nProxFacesMax * 3 * 3 * sizeof(double));
	cudaMalloc((void**)&comm->recvBuffGrad,   size * nProxFacesMax * 3 * 3 * sizeof(double));

	// Export struct
	cudaMalloc((void**)&exp->lcl2glbEl,    nElem  * sizeof(int));
	cudaMalloc((void**)&exp->lcl2glbNd,    nNodes * sizeof(int));
	cudaMalloc((void**)&exp->nElNds,       nElem  * sizeof(int));
	cudaMalloc((void**)&exp->elNd2lclNd,   nElem  * 8 * sizeof(int));
	cudaMalloc((void**)&exp->nEl4nd,       nNodes * sizeof(int));
	cudaMalloc((void**)&exp->idEl4nd,      nNodes * maxEl4nd * sizeof(int));
	cudaMalloc((void**)&exp->glbNdCord,    nNodes * 3 * sizeof(double));

	//

	cudaMemcpy(mesh->fc2el,       hLclData->fc2el,       nFcs * 2 * sizeof(int),    cudaMemcpyHostToDevice);
	free(hLclData->fc2el);

	cudaMemcpy(mesh->boundCond,   hLclData->boundCond,   nFcs * sizeof(int),        cudaMemcpyHostToDevice);
	free(hLclData->boundCond);

	cudaMemcpy(mesh->cellCenter,  hLclData->cellCenter,  nElem * 3 * sizeof(double), cudaMemcpyHostToDevice);
	free(hLclData->cellCenter);
	
	cudaMemcpy(mesh->fcCenter,    hLclData->fcCenter,    nFcs  * 3 * sizeof(double), cudaMemcpyHostToDevice);
	free(hLclData->fcCenter);

	cudaMemcpy(mesh->n,           hLclData->n,           nFcs  * 3 * sizeof(double), cudaMemcpyHostToDevice);
	free(hLclData->n);

	cudaMemcpy(mesh->nt1,         hLclData->nt1,         nFcs  * 3 * sizeof(double), cudaMemcpyHostToDevice);
	free(hLclData->nt1);

	cudaMemcpy(mesh->nt2,         hLclData->nt2,         nFcs  * 3 * sizeof(double), cudaMemcpyHostToDevice);
	free(hLclData->nt2);

	cudaMemcpy(mesh->area,        hLclData->area,        nFcs  * sizeof(double),     cudaMemcpyHostToDevice);
	free(hLclData->area);

	cudaMemcpy(mesh->volume,      hLclData->volume,      (nElem + nGhostElem) * sizeof(double),     cudaMemcpyHostToDevice);
	free(hLclData->volume);

	cudaMemcpy(comm->neigRank4fc, hLclData->neigRank4fc, nFcs * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->neigRank4fc);

	cudaMemcpy(comm->lclProx4fc, hLclData->lclProx4fc,
	           nFcs * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lclProx4fc);
	
	cudaMemcpy(comm->lclFc2idRcv, hLclData->lclFc2idRcv,
	           nFcs * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lclFc2idRcv);

	cudaMemcpy(exp->lcl2glbEl, hLclData->lcl2glbEl,
	           nElem * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lcl2glbEl);

	cudaMemcpy(exp->lcl2glbNd, hLclData->lcl2glbNd,
	           nNodes * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lcl2glbNd);

	cudaMemcpy(exp->nElNds, hLclData->nElNds,
		  nElem * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->nElNds);

	cudaMemcpy(exp->elNd2lclNd, hLclData->elNd2lclNd,
		   nElem * 8 * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->elNd2lclNd);

	cudaMemcpy(exp->nEl4nd, hLclData->nEl4nd,
	           nNodes * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->nEl4nd);

	cudaMemcpy(exp->idEl4nd, hLclData->idEl4nd,
	           nNodes * maxEl4nd * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->idEl4nd);

	cudaMemcpy(exp->glbNdCord, hLclData->glbNdCord,
	           nNodes * 3 * sizeof(double), cudaMemcpyHostToDevice);
	free(hLclData->glbNdCord);

	mesh->nElem         = nElem;
	mesh->nFcs          = nFcs;
	mesh->nOuter        = hLclData->nOuter;
	mesh->nProxTot      = hLclData->nProxTot;
	mesh->nGhostElem    = hLclData->nGhostElem;
	comm->nNeigRanks    = nNeigRanks;
	comm->nProxFacesMax = nProxFacesMax;
	comm->size          = size;
	comm->rank          = rank;
	exp->nGlbElem       = nGlbElem;
	exp->nGlbNodes      = nGlbNodes;
	exp->nNodes         = nNodes;
	exp->nElem          = nElem;
	exp->maxEl4nd       = maxEl4nd;
	
	comm->neigRanks  = (int*)malloc(nNeigRanks*sizeof(int));
	for(int i = 0; i < comm->nNeigRanks; i++){
		comm->neigRanks[i] = hLclData->neigRanks[i];
	}
	free(hLclData->neigRanks);
	
	comm->nProxFaces = (int*)malloc(nNeigRanks*sizeof(int));
	for(int i = 0; i < size; i++){
		comm->nProxFaces[i] = hLclData->nProxFaces[i];
	}
	free(hLclData->nProxFaces);

	comm->recvRequests = (MPI_Request**)malloc(3 * sizeof(MPI_Request *));
	for(int i = 0; i < 3; i++){
		comm->recvRequests[i] = (MPI_Request*)malloc(nNeigRanks * sizeof(MPI_Request));
	}
	
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("CUDA error in start device: %s\n", cudaGetErrorString(err));
	}
	
	free(hLclData);	 
	
}
	
