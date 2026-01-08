#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"
#include "bc.h"
#include "cons2prim.h"
#include "waveSpeeds.h"
#include "flux.h"
#include "starRegion.h"
#include "prim2cons.h"

__global__ void initializeRHS(int nElem, float *RHS){

	int i = threadIdx.x + blockIdx.x * blockDim.x;
	
	while(i < nElem*5){

		RHS[i] = 0.0;

		i += blockDim.x * gridDim.x;
	}
}

__global__ void computeFlux(INPUT *d_input, int start, int nFcs, int *fc2el, int *boundCond, float *n, float *nt1, float *nt2, float *area, float *primVecF, float *RHS){

	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;
	float solL[5];
	float solR[5];
	float vecS[5];
	float consF[5];
	int id = 0;

	int isLeftState=0, isStarRegion=0;
	int is_bc5,bc;

	float nVec[3],nt1Vec[3],nt2Vec[3];
	float areaF;
	float fluxAdv[15] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        float fluxHLLC[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};
	float uL_lcl[3], uR_lcl[3], c[3], cK;
	
	while(ifc < nFcs){	
		bc = boundCond[ifc];
		id = ( bc == 0 );
		
		solL[0] = primVecF[ifc*2*5 + 0*5 + 0];	
		solL[1] = primVecF[ifc*2*5 + 0*5 + 1];	
		solL[2] = primVecF[ifc*2*5 + 0*5 + 2];	
		solL[3] = primVecF[ifc*2*5 + 0*5 + 3];	
		solL[4] = primVecF[ifc*2*5 + 0*5 + 4];	
		
		solR[0] = primVecF[ifc*2*5 + id*5 + 0];	
		solR[1] = primVecF[ifc*2*5 + id*5 + 1];	
		solR[2] = primVecF[ifc*2*5 + id*5 + 2];	
		solR[3] = primVecF[ifc*2*5 + id*5 + 3];	
		solR[4] = primVecF[ifc*2*5 + id*5 + 4];	
				
		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];
		
		nt1Vec[0] = nt1[ifc*3 + 0];
		nt1Vec[1] = nt1[ifc*3 + 1];
		nt1Vec[2] = nt1[ifc*3 + 2];

		nt2Vec[0] = nt2[ifc*3 + 0];
		nt2Vec[1] = nt2[ifc*3 + 1];
		nt2Vec[2] = nt2[ifc*3 + 2];
		if(bc != 0){bcValues(d_input,nVec,nt1Vec,nt2Vec,solR,bc);}
		
		waveSpeeds(d_input,nVec,nt1Vec,nt2Vec,solL,solR,c,uL_lcl,uR_lcl);
		is_bc5 = (bc == 5);
		c[1] = (1-is_bc5)*c[1];
		
		isLeftState     = (c[1] >= 0 );
		isStarRegion    = (c[0] < 0 && c[2] >= 0);
		float (&sol)[5]    = isLeftState ? solL   : solR;
		float (&velLcl)[3] = isLeftState ? uL_lcl : uR_lcl;
		cK  = isLeftState ? c[0]   : c[2];

		if(isStarRegion){
			prim2cons(d_input,sol,consF);
			starRegion(d_input,sol,velLcl,cK,c[1],nVec,nt1Vec,nt2Vec,vecS); //define 2 arrays * 3 
			for(int iv = 0; iv < 5; iv++){
				fluxHLLC[iv] = cK * (vecS[iv]-consF[iv]);
			}
		}	
		flux(d_input,sol,fluxAdv,nVec);
			
		bc = boundCond[ifc];
		areaF = area[ifc];	
		for(int iv = 0; iv < 5; iv++){	
			atomicAdd( &RHS[fc2el[ifc*2+0]*5+iv], -(fluxAdv[iv*3 + 0]  + fluxAdv[iv*3 + 1] + fluxAdv[iv*3 + 2] + fluxHLLC[iv])*areaF);
			if(bc==0){atomicAdd( &RHS[fc2el[ifc*2+1]*5+iv], (fluxAdv[iv*3 + 0]  + fluxAdv[iv*3 + 1] + fluxAdv[iv*3 + 2] + fluxHLLC[iv])*areaF);}
		}
		
		ifc += blockDim.x * gridDim.x;
	}
}

__global__ void computeFluxProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclFc2idRcv, float *n, float *nt1, float *nt2, float *area, float *primVecF, float *recvBuff, float *RHS){

	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;
	float solL[5];
	float solR[5];
	float vecS[5];
	float consF[5];

	int isLeftState=0, isStarRegion=0;
	int is_bc5,bc;
	int neigRank, iProx;

	float nVec[3],nt1Vec[3],nt2Vec[3];
	float areaF;
	float fluxAdv[15] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        float fluxHLLC[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};
	float uL_lcl[3], uR_lcl[3], c[3], cK;

	while(ifc < nFcs){	
		
		solL[0] = primVecF[ifc*2*5 + 0*5 + 0];	
		solL[1] = primVecF[ifc*2*5 + 0*5 + 1];	
		solL[2] = primVecF[ifc*2*5 + 0*5 + 2];	
		solL[3] = primVecF[ifc*2*5 + 0*5 + 3];	
		solL[4] = primVecF[ifc*2*5 + 0*5 + 4];	
		
		neigRank = neigRank4fc[ifc];
		iProx    = lclFc2idRcv[ifc];
		solR[0] = recvBuff[neigRank*max*5 + iProx*5 + 0]; 
		solR[1] = recvBuff[neigRank*max*5 + iProx*5 + 1];
		solR[2] = recvBuff[neigRank*max*5 + iProx*5 + 2];
		solR[3] = recvBuff[neigRank*max*5 + iProx*5 + 3];
		solR[4] = recvBuff[neigRank*max*5 + iProx*5 + 4];
				
		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];
		
		nt1Vec[0] = nt1[ifc*3 + 0];
		nt1Vec[1] = nt1[ifc*3 + 1];
		nt1Vec[2] = nt1[ifc*3 + 2];

		nt2Vec[0] = nt2[ifc*3 + 0];
		nt2Vec[1] = nt2[ifc*3 + 1];
		nt2Vec[2] = nt2[ifc*3 + 2];
		
		waveSpeeds(d_input,nVec,nt1Vec,nt2Vec,solL,solR,c,uL_lcl,uR_lcl);
		is_bc5 = (bc == 5);
		c[1] = (1-is_bc5)*c[1];
		
		isLeftState     = (c[1] >= 0 );
		isStarRegion    = (c[0] < 0 && c[2] >= 0);
		float (&sol)[5]    = isLeftState ? solL   : solR;
		float (&velLcl)[3] = isLeftState ? uL_lcl : uR_lcl;
		cK  = isLeftState ? c[0]   : c[2];

		if(isStarRegion){
			prim2cons(d_input,sol,consF);
			starRegion(d_input,sol,velLcl,cK,c[1],nVec,nt1Vec,nt2Vec,vecS); 
			for(int iv = 0; iv < 5; iv++){
				fluxHLLC[iv] = cK * (vecS[iv]-consF[iv]);
			}
		}	
		flux(d_input,sol,fluxAdv,nVec);
			
		areaF = area[ifc];	
		for(int iv = 0; iv < 5; iv++){	
			atomicAdd( &RHS[fc2el[ifc*2+0]*5+iv], -(fluxAdv[iv*3 + 0]  + fluxAdv[iv*3 + 1] + fluxAdv[iv*3 + 2] + fluxHLLC[iv])*areaF);
		}
		
		ifc += blockDim.x * gridDim.x;
	}
}

void computeRHS(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evInitRHS, cudaEvent_t evFluxDone){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);

	initializeRHS<<<(mesh->nElem + 511)/512,512,0,sMain>>>(mesh->nElem, field->RHS);
	cudaEventRecord(evInitRHS, sMain);

	computeFlux<<<(nInner + 511)/512,512,0,sMain>>>(d_input,0, nInner, mesh->fc2el, mesh->boundCond, mesh->n, mesh->nt1, mesh->nt2, mesh->area, field->primVecF, field->RHS);
		
	computeViscousStresses<<<(nInner + 511)/512,512,0,sMain>>>(d_input, nInner, mesh->fc2el, mesh->boundCond, field->primVecF, field->grad, mesh->n, mesh->area, mesh->volume, field->RHS);
	
	MPI_Waitall(comm->nNeigRanks, comm->recvRequests[1], MPI_STATUSES_IGNORE);
	cudaStreamWaitEvent(sProx, evInitRHS, 0);
	cudaStreamSynchronize(sProx);   
	computeFluxProx<<<(nProxTh + 511)/512,512,0,sProx>>>(d_input,nInner, nFcs,comm->nProxFacesMax, mesh->fc2el, comm->neigRank4fc, comm->lclFc2idRcv, mesh->n, mesh->nt1, mesh->nt2, mesh->area, field->primVecF, comm->recvbuff, field->RHS);

	MPI_Waitall(comm->nNeigRanks, comm->recvRequests[2], MPI_STATUSES_IGNORE);
	cudaStreamSynchronize(sProx); 
	computeViscousStressesProx<<<(nProxTh + 511)/512,512,0,sProx>>>(d_input, nInner, nFcs, comm->nProxFacesMax, mesh->fc2el, comm->neigRank4fc, comm->lclFc2idRcv, field->primVecF, field->grad, mesh->n, mesh->area, mesh->volume, field->RHS, comm->recvBuffGrad, comm->recvbuff);
	cudaEventRecord(evFluxDone, sProx);

}
