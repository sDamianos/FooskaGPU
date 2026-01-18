#include "strdata.h"

__global__ void integrateConsVar(int iVar, int nElem, const double *volume, const double *consVec, double *blockSum){

	__shared__ double intVarBlock;
	
	if (threadIdx.x == 0){
		intVarBlock = 0.0;
	}

	__syncthreads();

	int iel = threadIdx.x + blockIdx.x * blockDim.x;
	double localVal = 0.0;

	double var, dVar; 
	while (iel < nElem) {
		if(iVar == 0){
			var  = 1.0 - consVec[iel * 10 + iVar];
		}else{
			var  = consVec[iel * 10 + iVar];
		}
		dVar = var * volume[iel];
		localVal += dVar;
		iel += blockDim.x * gridDim.x;
	}

	atomicAdd(&intVarBlock, localVal);
	__syncthreads();

	if (threadIdx.x == 0){
		blockSum[blockIdx.x] = intVarBlock;
	}
}

void consCalc(MESH *mesh, FIELD *field, cudaStream_t sMain){

	int nBlocks = (mesh->nElem + 511) / 512;
	int blockSize = 512;

	double *d_blockSum;
	cudaMalloc(&d_blockSum, nBlocks * sizeof(double));

	integrateConsVar<<<nBlocks, blockSize, 0, sMain>>>(input.consCalcIvar,mesh->nElem,mesh->volume,field->consVec,d_blockSum);

	cudaStreamSynchronize(sMain);

	double *h_blockSum = (double*) malloc(nBlocks * sizeof(double));
	cudaMemcpy(h_blockSum, d_blockSum, nBlocks * sizeof(double),cudaMemcpyDeviceToHost);
	
	cudaFree(d_blockSum);
	
	double intVarLcl = 0.0;
	for (int i = 0; i < nBlocks; i++){
		intVarLcl += h_blockSum[i];
	}
	free(h_blockSum);

	double intVar=0;
	MPI_Allreduce(&intVarLcl,&intVar,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0){	
		const char *mode;
		if(iStep == 0){
			mode = "w";
		}else{
			mode = "a";
		}

		FILE *fp = fopen("consVar.dat", mode);
		fprintf(fp, "%le %le\n", physicalTime, intVar);
		fclose(fp);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
