#include "strdata.h"
#include "eqos.h"
#include "rotateVector.h"
#include "bc.h"
#include "cons2prim.h"
#include "prim2cons.h"

__global__ void initializeTheta(int nElem, float *theta){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while(i < nElem*5){
		theta[i] = 1e9;
		i += blockDim.x * gridDim.x;
	}
}
__device__ float atomicMinFloat(float* addr, float value) { 
	int* addr_i = (int*)addr; 
	int old = *addr_i, assumed; 
	do { 
		assumed = old; 
		float old_val = __int_as_float(assumed); 
		if (old_val <= value) break; 
		old = atomicCAS(addr_i, assumed, __float_as_int(value)); 
	} while (assumed != old); 
	return __int_as_float(old); 
}

__global__ void computeLimiter(INPUT *d_input,int nFcs, int *fc2el, int *boundCond, float *primVecF, float *grad, float *cellCenter, float *faceCenter, float *wMin, float *wMax, float *extensionVec, float *theta){


	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/(2*5);
	int id  = (tid-ifc*5*2)/5;
	int iv  = (tid-ifc)%5; 

	float sol;
	float gradReg[3];
	float theta_temp;
	float xc,yc,zc,xf,yf,zf;

	float wmax, wmin, w_unlim;
	float phi;
	
	int bc;

	while(ifc < nFcs){
		
		bc       = boundCond[ifc];
		if(id != 1 || bc == 0){
		
			sol = primVecF[ifc*2*5 + id*5 + iv];
		
			for(int idr = 0; idr < 3; idr++){
				gradReg[idr] =  grad[fc2el[ifc*2+id]*5*3 + iv*3 + idr];
 			}

			xc = cellCenter[fc2el[ifc*2+id]*3 + 0];
			yc = cellCenter[fc2el[ifc*2+id]*3 + 1];
			zc = cellCenter[fc2el[ifc*2+id]*3 + 2];

			xf = faceCenter[ifc*3 + 0];
			yf = faceCenter[ifc*3 + 1];
			zf = faceCenter[ifc*3 + 2];
	
			w_unlim = sol + (xf-xc)*gradReg[0] + (yf-yc)*gradReg[1] + (zf-zc)*gradReg[2];
			
			extensionVec[ifc*5*2 + id*5 + iv] = (xf-xc)*gradReg[0] + (yf-yc)*gradReg[1] + (zf-zc)*gradReg[2];

			wmax = wMax[fc2el[ifc*2+id]*5 + iv]; 
			wmin = wMin[fc2el[ifc*2+id]*5 + iv];
	
			if(w_unlim-sol !=0 ){
				phi = (w_unlim > sol) ? (wmax - sol)/(2*(w_unlim-sol)) : (wmin - sol)/(2*(w_unlim-sol));
			}
			else{
				phi = 1;
			}
			theta_temp = fmaxf(0,  fminf(phi,1) ); //Minmode limieter
			
			atomicMinFloat(&theta[fc2el[ifc*2+id]*5 + iv],theta_temp);
		}
		tid += blockDim.x * gridDim.x;
	 	ifc = tid/(2*5);
		id  = (tid-ifc*5*2)/5;
		iv  = (tid-ifc)%5; 
	}
}


__global__ void updateFaces(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, float *theta, float *extensionVec, float *primVecF){
	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/(2*5);
	int id  = (tid-ifc*5*2)/5;
	int iv  = (tid-ifc)%5; 


	int bc;

	while(ifc < nFcs){
		bc = boundCond[ifc];
		if(id != 1 || bc == 0){	
			primVecF[ifc*2*5 + id*5 + iv] += theta[fc2el[ifc*2+id]*5 + iv]*extensionVec[ifc*5*2 + id*5 + iv];
		}
		tid += blockDim.x * gridDim.x;
	 	ifc = tid/(2*5);
		id  = (tid-ifc*5*2)/5;
		iv  = (tid-ifc)%5; 
	}
}

__global__ void updateFacesProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclProx4fc, float *theta, float *extensionVec, float *primVecF, float *sendBuff){
	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = start + tid/5;
	int iv  = tid%5; 

	int iProx,neigRank;

	while(ifc < nFcs){
		
		neigRank = neigRank4fc[ifc];
		iProx    = lclProx4fc[ifc];
					
		primVecF[ifc*2*5 + 0*5 + iv] += theta[fc2el[ifc*2+0]*5 + iv]*extensionVec[ifc*5*2 + 0*5 + iv];

		sendBuff[neigRank*max*5 + iProx*5 + iv] = primVecF[ifc*2*5 + 0*5 + iv];

		tid += blockDim.x * gridDim.x;
	 	ifc += tid/5; 
		iv   = tid%5; 
	}
}

void computeFaceValues(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evGreenGauss, cudaEvent_t evCompLimiter){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);
		
	initializeTheta<<<(mesh->nElem*5+511)/512,512,0,sMain>>>(mesh->nElem, field->theta);

	cudaStreamWaitEvent(sMain, evGreenGauss, 0);
	computeLimiter<<<(nFcs*2*5+511)/512,512,0,sMain>>>(d_input,nFcs,mesh->fc2el,mesh->boundCond,field->primVecF,field->grad,mesh->cellCenter,mesh->fcCenter,field->wMin,field->wMax,field->extensionVec,field->theta);
	cudaEventRecord(evCompLimiter,sMain);

	cudaStreamWaitEvent(sProx, evCompLimiter,0);
	updateFacesProx<<<(nProxTh*5+511)/512,512,0,sProx>>>(d_input,nInner,nFcs,comm->nProxFacesMax,mesh->fc2el,comm->neigRank4fc,comm->lclProx4fc,field->theta,field->extensionVec,field->primVecF, comm->sendbuff);
	
	updateFaces<<<(nInner*2*5+511)/512,512,0,sMain>>>(d_input,nInner,mesh->fc2el,mesh->boundCond,field->theta,field->extensionVec,field->primVecF);
	
	cudaStreamSynchronize(sProx);
	sendArray(comm,1,5,comm->sendbuff);
}

