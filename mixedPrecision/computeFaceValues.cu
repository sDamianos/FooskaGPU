#include "strdata.h"


__global__ void initializeTheta(int nElem, float *theta){

	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while(i < nElem*10){
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

__global__ void computeLimiter(INPUT *d_input,int nFcs, int *fc2el, int *boundCond, double *primVecF, float *grad, float *cellCenter, float *faceCenter, double *wMin, double *wMax, float *extensionVec, float *theta){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/(2*10);
	int id  = (tid-ifc*10*2)/10;
	int iv  = (tid-ifc)%10; 

	double sol;
	float gradReg[3];
	float theta_temp;
	float xc,yc,zc,xf,yf,zf;

	double w_unlim;
	float phi;
	double wmin,wmax;
	
	int bc;

	while(ifc < nFcs){
		
		bc       = boundCond[ifc];
		if(id != 1 || bc == 0){
		
			sol = primVecF[ifc*2*10 + id*10 + iv];
		
			for(int idr = 0; idr < 3; idr++){
				gradReg[idr] =  grad[fc2el[ifc*2+id]*10*3 + iv*3 + idr];
 			}
			//if(id==0 && iv == 6 && ifc == 68377){printf("iel %d ifc %d grad %e\n",fc2el[ifc*2+id], ifc, gradReg[0]);}

			xc = cellCenter[fc2el[ifc*2+id]*3 + 0];
			yc = cellCenter[fc2el[ifc*2+id]*3 + 1];
			zc = cellCenter[fc2el[ifc*2+id]*3 + 2];

			xf = faceCenter[ifc*3 + 0];
			yf = faceCenter[ifc*3 + 1];
			zf = faceCenter[ifc*3 + 2];
	
			w_unlim = sol + (xf-xc)*gradReg[0] + (yf-yc)*gradReg[1] + (zf-zc)*gradReg[2];
			
			extensionVec[ifc*10*2 + id*10 + iv] = (xf-xc)*gradReg[0] + (yf-yc)*gradReg[1] + (zf-zc)*gradReg[2];
			//if(id==0 && iv == 6 && ifc == 68377){printf("ifc %d  extension %e\n", ifc, extensionVec[ifc*10*2 + id*10 + iv]);}

			wmax = wMax[fc2el[ifc*2+id]*10 + iv]; 
			wmin = wMin[fc2el[ifc*2+id]*10 + iv];
//			if(fc2el[ifc*2+id]==17024 && iv == 6 && ifc == 68377){printf("ifc %d wmax wmin %e %e\n", ifc, wmax, wmin);}
	
			if(w_unlim-sol !=0 ){
				phi = (w_unlim > sol) ? (wmax - sol)/(2*(w_unlim-sol)) : (wmin - sol)/(2*(w_unlim-sol));
			}
			else{
				phi = 1;
			}
			theta_temp = fmax(0,  fmin(phi,1) ); //Minmode limieter
			//if(fc2el[ifc*2+id]==17023 && iv == 6){printf("ifc %d thetaF %le\n", ifc, theta_temp);}
			
			atomicMinFloat(&theta[fc2el[ifc*2+id]*10 + iv],theta_temp);
		}
		tid += blockDim.x * gridDim.x;
	 	ifc = tid/(2*10);
		id  = (tid-ifc*10*2)/10;
		iv  = (tid-ifc)%10; 
	}
}


__global__ void updateFaces(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, float *theta, float *extensionVec, double *primVecF){
	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int ifc = tid/(2*10);
	int id  = (tid-ifc*10*2)/10;
	int iv  = (tid-ifc)%10; 


	int bc;

	while(ifc < nFcs){
		bc = boundCond[ifc];
		
		if(id != 1 || bc == 0){	
			//if( ifc==68377 && iv==6 && id==0){printf("2nd Order: ifc %d theta %le extVec %le\n",ifc,theta[fc2el[ifc*2+id]*10 + iv],extensionVec[ifc*10*2 + id*10 + iv]);}
			primVecF[ifc*2*10 + id*10 + iv] += theta[fc2el[ifc*2+id]*10 + iv]*extensionVec[ifc*10*2 + id*10 + iv];
			//if( ifc==68377 && iv==6 && id==0){printf("2nd Order: ifc %d id %d primVec %.10le\n",ifc,id,primVecF[ifc*2*10 + id*10 + iv]);}
		}
		
		tid += blockDim.x * gridDim.x;
	 	ifc = tid/(2*10);
		id  = (tid-ifc*10*2)/10;
		iv  = (tid-ifc)%10; 
	}
}

__global__ void updateFacesProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclProx4fc, float *theta, float *extensionVec, double *primVecF, float *sendBuff){
	
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

//	cudaStreamWaitEvent(sMain, evGreenGauss, 0);
	computeLimiter<<<(nFcs*2*10+511)/512,512,0,sMain>>>(d_input,nFcs,mesh->fc2el,mesh->boundCond,field->primVecF,field->grad,mesh->cellCenter,mesh->fcCenter,field->wMin,field->wMax,field->extensionVec,field->theta);
//	cudaEventRecord(evCompLimiter,sMain);


//	cudaStreamWaitEvent(sProx, evCompLimiter,0);
//	updateFacesProx<<<(nProxTh*10+511)/512,512,0,sProx>>>(d_input,nInner,nFcs,comm->nProxFacesMax,mesh->fc2el,comm->neigRank4fc,comm->lclProx4fc,field->theta,field->extensionVec,field->primVecF, comm->sendbuff);
	
	updateFaces<<<(nInner*2*10+511)/512,512,0,sMain>>>(d_input,nInner,mesh->fc2el,mesh->boundCond,field->theta,field->extensionVec,field->primVecF);
	
//	cudaEventRecord(stop, sMain);
//	cudaEventSynchronize(stop);
 //	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("FaceValues time %3.5f ms\n", elapsedTime);
	
//	cudaStreamSynchronize(sProx);
//	sendArray(comm,1,10,comm->sendbuff);
}
