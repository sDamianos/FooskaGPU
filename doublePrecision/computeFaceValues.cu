#include "strdata.h"


__global__ void initializeTheta(int nElem, double *theta){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while(i < nElem*10){
		theta[i] = 1e9;
		i += blockDim.x * gridDim.x;
	}
}

__device__ double atomicMinDouble(double* addr, double value) {

	unsigned long long int* addr_ull = (unsigned long long int*)addr;
	unsigned long long int old = *addr_ull, assumed;

	do {
		assumed = old;
		double old_val = __longlong_as_double(assumed);
		if (old_val <= value) break;  // Already smaller, no update needed
		old = atomicCAS(addr_ull, assumed, __double_as_longlong(value));
	} while (assumed != old);
	
	return __longlong_as_double(old);  // returns previous value
}

__global__ void computeLimiter(INPUT *d_input,int nFcs, int *fc2el, int *boundCond, double *primVecF, double *grad, double *cellCenter, double *faceCenter, double *wMin, double *wMax, double *extensionVec, double *theta){


	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/(2*10);
	int id  = (tid-ifc*10*2)/10;
	int iv  = (tid-ifc)%10; 

	double sol;
	double gradReg[3];
	double theta_temp;
	double xc,yc,zc,xf,yf,zf;

	double wmax, wmin, w_unlim;
	double phi;
	
	int bc;

	while(ifc < nFcs){
		
		bc       = boundCond[ifc];
		if(id != 1 || bc == 0){
		
			sol = primVecF[ifc*2*10 + id*10 + iv];
		
			for(int idr = 0; idr < 3; idr++){
				gradReg[idr] =  grad[fc2el[ifc*2+id]*10*3 + iv*3 + idr];
 			}

			xc = cellCenter[fc2el[ifc*2+id]*3 + 0];
			yc = cellCenter[fc2el[ifc*2+id]*3 + 1];
			zc = cellCenter[fc2el[ifc*2+id]*3 + 2];

			xf = faceCenter[ifc*3 + 0];
			yf = faceCenter[ifc*3 + 1];
			zf = faceCenter[ifc*3 + 2];
	
			w_unlim = sol + (xf-xc)*gradReg[0] + (yf-yc)*gradReg[1] + (zf-zc)*gradReg[2];
			
			extensionVec[ifc*10*2 + id*10 + iv] = (xf-xc)*gradReg[0] + (yf-yc)*gradReg[1] + (zf-zc)*gradReg[2];

			wmax = wMax[fc2el[ifc*2+id]*10 + iv]; 
			wmin = wMin[fc2el[ifc*2+id]*10 + iv];
	
			if(w_unlim-sol !=0 ){
				phi = (w_unlim > sol) ? (wmax - sol)/(2*(w_unlim-sol)) : (wmin - sol)/(2*(w_unlim-sol));
			}
			else{
				phi = 1;
			}
			theta_temp = fmax(0,  fmin(phi,1) ); //Minmode limieter
			
			atomicMinDouble(&theta[fc2el[ifc*2+id]*10 + iv],theta_temp);
		}
		tid += blockDim.x * gridDim.x;
	 	ifc = tid/(2*10);
		id  = (tid-ifc*10*2)/10;
		iv  = (tid-ifc)%10; 
	}
}


__global__ void updateFaces(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, double *theta, double *extensionVec, double *primVecF){
	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/(2*10);
	int id  = (tid-ifc*10*2)/10;
	int iv  = (tid-ifc)%10; 


	int bc;

	while(ifc < nFcs){
		bc = boundCond[ifc];
		
		if(id != 1 || bc == 0){	
			primVecF[ifc*2*10 + id*10 + iv] += theta[fc2el[ifc*2+id]*10 + iv]*extensionVec[ifc*10*2 + id*10 + iv];
		}
		
		tid += blockDim.x * gridDim.x;
	 	ifc = tid/(2*10);
		id  = (tid-ifc*10*2)/10;
		iv  = (tid-ifc)%10; 
	}
}

__global__ void updateFacesProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclProx4fc, double *theta, double *extensionVec, double *primVecF, double *sendBuff){
	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = start + tid/10;
	int iv  = tid%10; 

	int iProx,neigRank;

	while(ifc < nFcs){
		
		neigRank = neigRank4fc[ifc];
		iProx    = lclProx4fc[ifc];
					
		primVecF[ifc*2*10 + 0*10 + iv] += theta[fc2el[ifc*2+0]*10 + iv]*extensionVec[ifc*10*2 + 0*10 + iv];

		sendBuff[neigRank*max*10 + iProx*10 + iv] = primVecF[ifc*2*10 + 0*10 + iv];

		tid += blockDim.x * gridDim.x;
	 	ifc += tid/10; 
		iv   = tid%10; 
	}
}

void computeFaceValues(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evGreenGauss, cudaEvent_t evCompLimiter){
	
	int nInner  = mesh->nFcs-mesh->nProxTot;
	int nProx   = mesh->nProxTot;
	int nFcs    = mesh->nFcs;
	int nProxTh = max(1,nProx);

//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);	
//	cudaEventRecord(start,sMain);	
		
	initializeTheta<<<(mesh->nElem*10+511)/512,512,0,sMain>>>(mesh->nElem, field->theta);

	cudaStreamWaitEvent(sMain, evGreenGauss, 0);
	computeLimiter<<<(nFcs*2*10+511)/512,512,0,sMain>>>(d_input,nFcs,mesh->fc2el,mesh->boundCond,field->primVecF,field->grad,mesh->cellCenter,mesh->fcCenter,field->wMin,field->wMax,field->extensionVec,field->theta);
	cudaEventRecord(evCompLimiter,sMain);


	cudaStreamWaitEvent(sProx, evCompLimiter,0);
	updateFacesProx<<<(nProxTh*10+511)/512,512,0,sProx>>>(d_input,nInner,nFcs,comm->nProxFacesMax,mesh->fc2el,comm->neigRank4fc,comm->lclProx4fc,field->theta,field->extensionVec,field->primVecF, comm->sendbuff);
	
	updateFaces<<<(nInner*2*10+511)/512,512,0,sMain>>>(d_input,nInner,mesh->fc2el,mesh->boundCond,field->theta,field->extensionVec,field->primVecF);
	
//	cudaEventRecord(stop, sMain);
//	cudaEventSynchronize(stop);
 //	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("FaceValues time %3.5f ms\n", elapsedTime);
	
	cudaStreamSynchronize(sProx);
	sendArray(comm,1,10,comm->sendbuff);
}
