#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"
#include "bc.h"
#include "cons2prim.h"
#include "waveSpeeds.h"
#include "flux.h"
#include "starRegion.h"
#include "prim2cons.h"


__global__ void initializeRHS(int nElem, double *RHS){

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	
	while(i < nElem*10){

		RHS[i] = 0.0;

		i += blockDim.x * gridDim.x;
	}
}

__device__ void computeFluxHllc(double cK, double (&vecS)[10], double (&consF)[10], double (&fluxHLLC)[10]){
	
	fluxHLLC[0] = 0;
	fluxHLLC[1] = cK * (vecS[1]-consF[1]);
	fluxHLLC[2] = cK * (vecS[2]-consF[2]);
	fluxHLLC[3] = cK * (vecS[3]-consF[3]);
	fluxHLLC[4] = cK * (vecS[4]-consF[4]);
	fluxHLLC[5] = cK * (vecS[5]-consF[5]);
	fluxHLLC[6] = 0;
	fluxHLLC[7] = 0;
	fluxHLLC[8] = cK * (vecS[8]-consF[8]);
	fluxHLLC[9] = cK * (vecS[9]-consF[9]);
}

__device__ __noinline__ double computeCenterPress(INPUT *d_input, double *consVec){
	
	double avf;
	double rhoL,rhoG,eL,eG,pressL,pressG;

	avf    = fmin(fmax(d_input->avfMin/100,consVec[0]),1.0-d_input->avfMin/100);
	rhoL = consVec[1]/avf;
	rhoG = consVec[2]/(1-avf);
	eL  = consVec[6]/(avf*rhoL);
	eG  = consVec[7]/((1-avf)*rhoG); 
	pressL = eqos(d_input,0,0,rhoL,eL);
	pressG = eqos(d_input,1,0,rhoG,eG);

	return avf*pressL + (1-avf)*pressG;
}

__device__ __noinline__ double pressureCorrection(INPUT *d_input, double *sol, double *consVec){

	double pC = computeCenterPress(d_input,consVec);	 
	double pF = sol[0]*sol[6] + (1-sol[0])*sol[7];
	return pF - pC; 
}

__device__ __noinline__ void atomicAddRHS(int id, double *fluxAdv, double *fluxHLLC, double *fluxNonCons, double pressTerm, double *nVec, double areaF, double *RHS){
		
		double sign = (id == 0) ? -1.0 : 1.0;
			
		atomicAdd( &RHS[0], sign*(fluxAdv[0] - fluxNonCons[0*2+id])*areaF);
		atomicAdd( &RHS[1], sign*(fluxAdv[1] + fluxHLLC[1])*areaF);
		atomicAdd( &RHS[2], sign*(fluxAdv[2] + fluxHLLC[2])*areaF);
		atomicAdd( &RHS[3], sign*(fluxAdv[3] + fluxHLLC[3] + pressTerm*nVec[0])*areaF);
		atomicAdd( &RHS[4], sign*(fluxAdv[4] + fluxHLLC[4] + pressTerm*nVec[1])*areaF);
		atomicAdd( &RHS[5], sign*(fluxAdv[5] + fluxHLLC[5] + pressTerm*nVec[2])*areaF);
		atomicAdd( &RHS[6], sign*(fluxAdv[6] + fluxNonCons[1*2+id])*areaF);
		atomicAdd( &RHS[7], sign*(fluxAdv[7] + fluxNonCons[2*2+id])*areaF);
		atomicAdd( &RHS[8], sign*(fluxAdv[8] + fluxHLLC[8])*areaF);
		atomicAdd( &RHS[9], sign*(fluxAdv[9] + fluxHLLC[9])*areaF);
}
		
		
__global__ void computeFlux(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, double *n, double *nt1, double *nt2, double *area, double *primVecF, double *consVec, double *RHS){

	int ifc = threadIdx.x + blockIdx.x * blockDim.x;
	double solL[10];
	double solR[10];
	double vecS[10];
	double consF[10];
	int id = 0;

	int isLeftState=0, isStarRegion=0;
	int bc;

	double nVec[3],nt1Vec[3],nt2Vec[3];
	double areaF;
	double uNonCons[3];
	double fluxNonCons[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double fluxHLLC[10]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
	double fluxAdv[15]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
	double uL_lcl[3], uR_lcl[3], c[3], cK;
	
	while(ifc < nFcs){	
		bc = boundCond[ifc];
		id = ( bc == 0 );
		
		//if(fc2el[ifc*2+0]==3150538||fc2el[ifc*2+1]==3150538){printf("ifc %d ny*A %e\n",ifc, area[ifc]* n[ifc*3 +1]);}
		for(int iVar = 0; iVar < 10; iVar++){
			solL[iVar] = primVecF[ifc*2*10 + 0 *10 + iVar];	
			solR[iVar] = primVecF[ifc*2*10 + id*10 + iVar];	
		}
				
		for(int idr = 0; idr < 3; idr++){
			nVec[idr] = n[ifc*3 + idr];
			nt1Vec[idr] = nt1[ifc*3 + idr];
			nt2Vec[idr] = nt2[ifc*3 + idr];
		}	
		if(bc != 0){bcValues(d_input,nVec,nt1Vec,nt2Vec,solR,bc);}
		
		waveSpeeds(d_input,nVec,nt1Vec,nt2Vec,solL,solR,c,uL_lcl,uR_lcl);
		
		isLeftState        = (c[1] >= 0 );
		isStarRegion       = (c[0] < 0 && c[2] >= 0);
		double (&sol)[10]   = isLeftState ? solL   : solR;
		double (&velLcl)[3] = isLeftState ? uL_lcl : uR_lcl;
		cK  = isLeftState  ? c[0]   : c[2];

		uNonCons[0] = sol[3]; uNonCons[1] = sol[4]; uNonCons[2] = sol[5];

		prim2cons(d_input,sol,consF);
		if(isStarRegion){
			starRegion(d_input,sol,velLcl,cK,c[1],nVec,nt1Vec,nt2Vec,uNonCons,vecS,consF); //define 2 arrays * 3 
			computeFluxHllc(cK,vecS,consF,fluxHLLC);
		}	
		flux(d_input,isStarRegion,sol,fluxAdv,nVec,uNonCons,consF[6],consF[7]);
		fluxSource(d_input,uNonCons,&consVec[fc2el[ifc*2+0]*10],&consVec[fc2el[ifc*2+id]*10],nVec,fluxNonCons);
		
		double pressureCor[2];
		pressureCor[0] = pressureCorrection(d_input,sol,&consVec[fc2el[ifc*2+0]*10]);
		if(bc==0){pressureCor[1] = pressureCorrection(d_input,sol,&consVec[fc2el[ifc*2+1]*10]);}
		
		//if(fc2el[ifc*2+0]==3150538||fc2el[ifc*2+1]==3150538){printf("fluxAdv %e fluxHllc of V %e\n", fluxAdv[4], fluxHLLC[4]);}
			
		bc = boundCond[ifc];
		areaF = area[ifc];	
		atomicAddRHS(0,fluxAdv,fluxHLLC,fluxNonCons,pressureCor[0],nVec,areaF,&RHS[fc2el[ifc*2+0]*10]);
		if(bc==0){atomicAddRHS(1,fluxAdv,fluxHLLC,fluxNonCons,pressureCor[1],nVec,areaF,&RHS[fc2el[ifc*2+1]*10]);}
		
		ifc += blockDim.x * gridDim.x;
	}
}



__global__ void computeFluxProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclFc2idRcv, double *n, double *nt1, double *nt2, double *area, double *primVecF, double *recvBuff, double *consVec, double *RHS){

	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;
	double solL[10];
	double solR[10];
	double vecS[10];
	double consF[10];

	int isLeftState=0, isStarRegion=0;
	int neigRank, iProx;

	double nVec[3],nt1Vec[3],nt2Vec[3];
	double areaF;
	double uNonCons[3];
	double fluxNonCons[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double fluxHLLC[10]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
	double fluxAdv[15]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
	double uL_lcl[3], uR_lcl[3], c[3], cK;

	while(ifc < nFcs){	
		
		neigRank = neigRank4fc[ifc];
		iProx    = lclFc2idRcv[ifc];
		for(int iVar = 0; iVar < 10; iVar++){
			solL[iVar] = primVecF[ifc*2*10 + 0 *10 + iVar];	
			solR[iVar] = recvBuff[neigRank*max*10 + iProx*10 + iVar];
		}
		
		for(int idr = 0; idr < 3; idr++){
			nVec[idr] = n[ifc*3 + idr];
			nt1Vec[idr] = nt1[ifc*3 + idr];
			nt2Vec[idr] = nt2[ifc*3 + idr];
		}	
					
		waveSpeeds(d_input,nVec,nt1Vec,nt2Vec,solL,solR,c,uL_lcl,uR_lcl);
		
		isLeftState     = (c[1] >= 0 );
		isStarRegion    = (c[0] < 0 && c[2] >= 0);
		double (&sol)[10]   = isLeftState ? solL   : solR;
		double (&velLcl)[3] = isLeftState ? uL_lcl : uR_lcl;
		cK  = isLeftState ? c[0]   : c[2];
		
		uNonCons[0] = sol[3]; uNonCons[1] = sol[4]; uNonCons[2] = sol[5];

		prim2cons(d_input,sol,consF);	
		if(isStarRegion){
			starRegion(d_input,sol,velLcl,cK,c[1],nVec,nt1Vec,nt2Vec,uNonCons,vecS,consF); 
			computeFluxHllc(cK,vecS,consF,fluxHLLC);
		}	
		flux(d_input,isStarRegion,sol,fluxAdv,nVec,uNonCons,consF[6],consF[7]);
		fluxSource(d_input,uNonCons,&consVec[fc2el[ifc*2+0]*10],&consVec[fc2el[ifc*2+1]*10],nVec,fluxNonCons);

		double pressureCor[2];
		pressureCor[0] = pressureCorrection(d_input,sol,&consVec[fc2el[ifc*2+0]*10]);
				
		areaF = area[ifc];	
		atomicAddRHS(0,fluxAdv,fluxHLLC,fluxNonCons,pressureCor[0],nVec,areaF,&RHS[fc2el[ifc*2+0]*10]);

		ifc += blockDim.x * gridDim.x;
	}
	
}


void computeRHS(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitRHS, cudaEvent_t evFluxDone){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);
	
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);	
//	cudaEventRecord(start,sMain);	

	initializeRHS<<<(mesh->nElem*10 + 511)/512,512,0,sMain>>>(mesh->nElem, field->RHS);
	cudaEventRecord(evInitRHS, sMain);

	computeFlux<<<(nInner + 511)/512,512,0,sMain>>>(d_input, nInner, mesh->fc2el, mesh->boundCond, mesh->n, mesh->nt1, mesh->nt2, mesh->area, field->primVecF, field->consVec, field->RHS);
	
//	cudaEventRecord(stop, sMain);
//	cudaEventSynchronize(stop);
//	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("RHS time %3.5f ms\n", elapsedTime);
	
	if(input.viscFlag==1){	
		computeViscousStresses<<<(nInner + 511)/512,512,0,sMain>>>(d_input, nInner, mesh->fc2el, mesh->boundCond, field->primVecF, field->grad, mesh->n, mesh->area, mesh->volume, field->RHS);
	}
	
	MPI_Waitall(comm->nNeigRanks, comm->recvRequests[1], MPI_STATUSES_IGNORE);
	cudaStreamWaitEvent(sProx, evInitRHS, 0);
	cudaStreamSynchronize(sProx);   
	computeFluxProx<<<(nProxTh + 511)/512,512,0,sProx>>>(d_input,nInner, nFcs,comm->nProxFacesMax, mesh->fc2el, comm->neigRank4fc, comm->lclFc2idRcv, mesh->n, mesh->nt1, mesh->nt2, mesh->area, field->primVecF, comm->recvbuff, field->consVec, field->RHS);

	if(input.viscFlag==1){
		MPI_Waitall(comm->nNeigRanks, comm->recvRequests[2], MPI_STATUSES_IGNORE);
		cudaStreamSynchronize(sProx); 
		computeViscousStressesProx<<<(nProxTh + 511)/512,512,0,sProx>>>(d_input, nInner, nFcs, comm->nProxFacesMax, mesh->fc2el, comm->neigRank4fc, comm->lclFc2idRcv, field->primVecF, field->grad, mesh->n, mesh->area, mesh->volume, field->RHS, comm->recvBuffGrad, comm->recvbuff);
	}
	cudaEventRecord(evFluxDone, sProx);

}
