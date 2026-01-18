#include "strdata.h"

__global__ void rk3stage0(int nElem, double *consVec0, double *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/10;
	int iVar  = tid%10;
	
	while(iel < nElem){
		
		consVec0[iel*10 + iVar] = consVec[iel*10 + iVar];
		
		tid += blockDim.x * gridDim.x;
		iel = tid/10;
		iVar  = tid%10;
	}
}

__global__ void rk3stage1(INPUT *d_input, int nElem, double dt, double *volume, double *RHS, double *consVec0, double *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/10;
	int iVar  = tid%10;
	int neglect;

	neglect = (iVar == 3 && d_input->neglectX == 1); 
	neglect = (iVar == 4 && d_input->neglectY == 1); 
	neglect = (iVar == 5 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*10 + iVar] = (1-neglect) * ( consVec0[iel*10 + iVar] + dt * (1.0/volume[iel]) * RHS[iel*10+iVar] );
		
		tid += blockDim.x * gridDim.x;
		iel = tid/10;
		iVar  = tid%10;
	}
}

__global__ void rk3stage2(INPUT *d_input, int nElem, double dt, double *volume, double *RHS, double *consVec0, double *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/10;
	int iVar  = tid%10;
	int neglect;

	neglect = (iVar == 3 && d_input->neglectX == 1); 
	neglect = (iVar == 4 && d_input->neglectY == 1); 
	neglect = (iVar == 5 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*10 + iVar] = (1-neglect) * ( (3.0/4.0)*consVec0[iel*10 + iVar] + (1.0/4.0)*consVec[iel*10 + iVar] + (1.0/4.0) * dt * (1.0/volume[iel]) * RHS[iel*10+iVar] );
		
		tid += blockDim.x * gridDim.x;
		iel = tid/10;
		iVar  = tid%10;
	}
}

__global__ void rk3stage3(INPUT *d_input, int nElem, double dt, double *volume, double *RHS, double *consVec0, double *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/10;
	int iVar  = tid%10;
	int neglect;

	neglect = (iVar == 3 && d_input->neglectX == 1); 
	neglect = (iVar == 4 && d_input->neglectY == 1); 
	neglect = (iVar == 5 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*10 + iVar] = (1-neglect) * ( (1.0/3.0)*consVec0[iel*10 + iVar] + (2.0/3.0)*consVec[iel*10 + iVar] +  (2.0/3.0) * dt * (1.0/volume[iel]) * RHS[iel*10+iVar] );
		
		tid += blockDim.x * gridDim.x;
		iel = tid/10;
		iVar  = tid%10;
	}
}


__global__ void rungeKutta(INPUT *d_input, int nElem, double dt, double *volume, double *RHS, double *consVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/10;
	int iVar  = tid%10;
	int neglect;

	neglect = (iVar == 3 && d_input->neglectX == 1); 
	neglect = (iVar == 4 && d_input->neglectY == 1); 
	neglect = (iVar == 5 && d_input->neglectZ == 1);

	while(iel < nElem){
		
		consVec[iel*10 + iVar] = (1-neglect) * ( consVec[iel*10 + iVar] + dt * (1.0/volume[iel]) * RHS[iel*10+iVar] );
		//if(iel==3150538&&iVar==4){printf("RHS of var %d is %e\n", iVar, dt * (1.0/volume[iel]) * RHS[iel*10+iVar]);}
		
		tid += blockDim.x * gridDim.x;
		iel = tid/10;
		iVar  = tid%10;
	}
}

void update(int irk, MESH *mesh, FIELD *field, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evFluxDone, cudaEvent_t evStepDone){

	cudaStreamWaitEvent(sMain, evFluxDone, 0);
	
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);	
//	cudaEventRecord(start,sMain);	
	
	switch(input.rkSteps){
		case 1:
			rungeKutta<<<(mesh->nElem*10+511)/512,512,0,sMain>>>(d_input,mesh->nElem,input.dt,mesh->volume,field->RHS,field->consVec);
		break;

		case 3:
			switch(irk){
				case 0: rk3stage1<<<(mesh->nElem*10+511)/512,512,0,sMain>>>(d_input,mesh->nElem,input.dt,mesh->volume,field->RHS,field->consVec0,field->consVec); break;
				case 1: rk3stage2<<<(mesh->nElem*10+511)/512,512,0,sMain>>>(d_input,mesh->nElem,input.dt,mesh->volume,field->RHS,field->consVec0,field->consVec); break;
				case 2: rk3stage3<<<(mesh->nElem*10+511)/512,512,0,sMain>>>(d_input,mesh->nElem,input.dt,mesh->volume,field->RHS,field->consVec0,field->consVec); break;
			}
		break;
	}
//	cudaEventRecord(stop, sMain);
//	cudaEventSynchronize(stop);
  //  	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("update time %3.5f ms\n", elapsedTime);

	if(irk < input.rkSteps - 1){
		cudaEventRecord(evStepDone, sMain);
		cudaStreamWaitEvent(sProx, evStepDone, 0);
	}
	
}

		
void initRk(MESH *mesh, FIELD *field, cudaStream_t sMain){

		rk3stage0<<<(mesh->nElem*10+511)/512,512,0,sMain>>>(mesh->nElem, field->consVec0, field->consVec);

}
