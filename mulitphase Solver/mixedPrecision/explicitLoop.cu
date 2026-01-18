#include "strdata.h"

void explicitLoop(MESH *mesh, EXPORT *exp, COMM *comm, FIELD *field){
	
	int iExp = 0;

	cudaDeviceSynchronize();
	MPI_Barrier(MPI_COMM_WORLD);
			
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);	
	cudaEventRecord(start,0);	

	cudaEvent_t evInitGrads, evGreenGauss, evCompLimiter, evInitRHS, evFluxDone, evStepDone;
	cudaEventCreate(&evInitGrads); 
	cudaEventCreate(&evGreenGauss);
	cudaEventCreate(&evCompLimiter);
	cudaEventCreate(&evInitRHS);
	cudaEventCreate(&evFluxDone); 
	cudaEventCreate(&evStepDone);

	cudaStream_t sMain, sProx;
	cudaStreamCreateWithFlags(&sMain, cudaStreamNonBlocking);
	cudaStreamCreateWithFlags(&sProx, cudaStreamNonBlocking);

	iStep = 0;
	while(physicalTime < input.finalTime){
		physicalTime += input.dt;
		
		if(input.rkSteps>1){initRk(mesh, field, sMain);}
		for(int irk = 0; irk < input.rkSteps; irk++){
			recvArray(comm,0,10,comm->recvbuff);
			recvArray(comm,1,10,comm->recvbuff);
			recvArrayFloat(comm,2,3*3,comm->recvBuffGrad);	
		
			computePrimVec(mesh, field, comm, sMain, sProx);

			computeGrad(mesh, field, comm, sMain, sProx, evInitGrads, evGreenGauss);

			computeFaceValues(mesh, field, comm, sMain, sProx, evGreenGauss, evCompLimiter);

			computeRHS(mesh, field, comm, sMain, sProx, evInitRHS, evFluxDone);	

			if(irk==0&&iStep%input.adaptDtItr==0){computeTimeStep(mesh, field, comm, sMain, sProx);}

			update(irk, mesh, field, sMain, sProx, evFluxDone, evStepDone);
		}
		
		relaxation(mesh, field, comm,  sMain, sProx, evStepDone);
		
		if(input.consCalcFlag==1&&iStep%input.adaptDtItr==0){consCalc(mesh, field, sMain);}
		 	
		if(physicalTime >= input.exportTime * (iExp +1)){export_plt(exp,field); iExp++;}
		if((iStep+1)%input.restartStep==0){writeFieldRestart(comm->rank, mesh, field);}
		if((iStep+1)%1000 == 0 && comm->rank==0){
			cudaEventRecord(stop, 0);
			cudaEventSynchronize(stop);
    			float   elapsedTime;
			cudaEventElapsedTime(&elapsedTime,start,stop);
			printf("step %d Done Time %3.5f ms\n", iStep+1, elapsedTime);
		}
		iStep++;
	}
		
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("CUDA error: %s\n", cudaGetErrorString(err));
	}
}

