#include "strdata.h"

__global__ void rk3stage0(int nElem, float *consVec0, float *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/5;
	int iVar  = tid%5;
	
	while(iel < nElem){
		
		consVec0[iel*5 + iVar] = consVec[iel*5 + iVar];
		
		tid += blockDim.x * gridDim.x;
		iel = tid/5;
		iVar  = tid%5;
	}
}

__global__ void rk3stage1(INPUT *d_input, int nElem, float *volume, float *RHS, float *consVec0, float *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/5;
	int iVar  = tid%5;
	int neglect;

	neglect = (iVar == 1 && d_input->neglectX == 1); 
	neglect = (iVar == 2 && d_input->neglectY == 1); 
	neglect = (iVar == 3 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*5 + iVar] = (1-neglect) * ( consVec0[iel*5 + iVar] + d_input->dt * (1.0/volume[iel]) * RHS[iel*5+iVar] );
		
		tid += blockDim.x * gridDim.x;
		iel = tid/5;
		iVar  = tid%5;
	}
}

__global__ void rk3stage2(INPUT *d_input, int nElem, float *volume, float *RHS, float *consVec0, float *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/5;
	int iVar  = tid%5;
	int neglect;

	neglect = (iVar == 1 && d_input->neglectX == 1); 
	neglect = (iVar == 2 && d_input->neglectY == 1); 
	neglect = (iVar == 3 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*5 + iVar] = (1-neglect) * ( (3.0/4.0)*consVec0[iel*5 + iVar] + (1.0/4.0)*consVec[iel*5 + iVar] + (1.0/4.0) * d_input->dt * (1.0/volume[iel]) * RHS[iel*5+iVar] );
		
		tid += blockDim.x * gridDim.x;
		iel = tid/5;
		iVar  = tid%5;
	}
}

__global__ void rk3stage3(INPUT *d_input, int nElem, float *volume, float *RHS, float *consVec0, float *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/5;
	int iVar  = tid%5;
	int neglect;

	neglect = (iVar == 1 && d_input->neglectX == 1); 
	neglect = (iVar == 2 && d_input->neglectY == 1); 
	neglect = (iVar == 3 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*5 + iVar] = (1-neglect) * ( (1.0/3.0)*consVec0[iel*5 + iVar] + (2.0/3.0)*consVec[iel*5 + iVar] +  (2.0/3.0) * d_input->dt * (1.0/volume[iel]) * RHS[iel*5+iVar] );
		
		tid += blockDim.x * gridDim.x;
		iel = tid/5;
		iVar  = tid%5;
	}
}


__global__ void rungeKutta(INPUT *d_input, int nElem, float *volume, float *RHS, float *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/5;
	int iVar  = tid%5;
	int neglect;

	neglect = (iVar == 1 && d_input->neglectX == 1); 
	neglect = (iVar == 2 && d_input->neglectY == 1); 
	neglect = (iVar == 3 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*5 + iVar] = (1-neglect) * ( consVec[iel*5 + iVar] + d_input->dt * (1.0/volume[iel]) * RHS[iel*5+iVar] );
		
		tid += blockDim.x * gridDim.x;
		iel = tid/5;
		iVar  = tid%5;
	}
}

void update(int irk, MESH *mesh, FIELD *field, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evFluxDone, cudaEvent_t evStepDone){

	cudaStreamWaitEvent(sMain, evFluxDone, 0);
	
	switch(input.rkSteps){
		case 1:
			rungeKutta<<<(mesh->nElem*5+511)/512,512,0,sMain>>>(d_input,mesh->nElem,mesh->volume,field->RHS,field->consVec);
		break;

		case 3:
			switch(irk){
				case 0: rk3stage1<<<(mesh->nElem*5+511)/512,512,0,sMain>>>(d_input,mesh->nElem,mesh->volume,field->RHS,field->consVec0,field->consVec); break;
				case 1: rk3stage2<<<(mesh->nElem*5+511)/512,512,0,sMain>>>(d_input,mesh->nElem,mesh->volume,field->RHS,field->consVec0,field->consVec); break;
				case 2: rk3stage3<<<(mesh->nElem*5+511)/512,512,0,sMain>>>(d_input,mesh->nElem,mesh->volume,field->RHS,field->consVec0,field->consVec); break;
			}
		break;
	}
	cudaEventRecord(evStepDone, sMain);
	cudaStreamWaitEvent(sProx, evStepDone, 0);
}

		
void initRk(MESH *mesh, FIELD *field, cudaStream_t sMain){

		rk3stage0<<<(mesh->nElem*5+511)/512,512,0,sMain>>>(mesh->nElem, field->consVec0, field->consVec);

}
