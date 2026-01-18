#include "strdata.h"
#include "eqos.h"
//====== Phase 0 = LIQUID ======================== Phase 1 = GAS ============
__constant__ double press1 = 11.4e5;	 __constant__ double press2 = 2.6e5;
__constant__ double temp1  = 115.3;	 __constant__ double temp2  = 115.3;
__constant__ double u1     = 0;		 __constant__ double u2     = 0.0;
__constant__ double v1     = 0.0;	 __constant__ double v2     = 0.0;
__constant__ double w1     = 0.0;	 __constant__ double w2     = 0.0;
	
__global__ void initConsVec(int nElem, double *consVec){
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while(i<nElem*10){
		consVec[i]=1.0;
		i += gridDim.x * blockDim.x;
	}
}

__global__ void initVec(INPUT *d_input, int nElem, int NEQ, double *cellCenter, double *consVec){
	
	int iel = threadIdx.x + blockIdx.x * blockDim.x;
	double avfL,avfG,rhoL,rhoG,rhoMix,u,v,w,press,temp,eIntL,eIntG;
	int is_region_1 = 0;
	
	while(iel < nElem){
		
		is_region_1 = (cellCenter[iel*3 + 0] < 0.0 );
	
		avfL  = is_region_1 ? 1.0-d_input->avfMin : d_input->avfMin;
		press = is_region_1 ? press1 	          : press2;
		temp  = is_region_1 ? temp1  	          : temp2;
		u     = is_region_1 ? u1     	          : u2;
		v     = is_region_1 ? v1      	          : v2;
		w     = is_region_1 ? w1     	          : w2;
		
		avfG  = 1 - avfL; 
		
		rhoL  = eqos(d_input,0,2,press,temp);
		eIntL = eqos(d_input,0,3,press,rhoL);
		
		rhoG  = eqos(d_input,1,2,press,temp);
		eIntG  = eqos(d_input,1,3,press,rhoG);

		rhoMix = avfL*rhoL + avfG*rhoG;
		
		consVec[iel*NEQ + 0] = avfL;
		consVec[iel*NEQ + 1] = avfL*rhoL;
		consVec[iel*NEQ + 2] = avfG*rhoG;
		consVec[iel*NEQ + 3] = rhoMix*u;
		consVec[iel*NEQ + 4] = rhoMix*v;
		consVec[iel*NEQ + 5] = rhoMix*w;
		consVec[iel*NEQ + 6] = avfL*rhoL*eIntL;
		consVec[iel*NEQ + 7] = avfG*rhoG*eIntG;
		consVec[iel*NEQ + 8] = avfL*rhoL*eIntL + avfG*rhoG*eIntG + 0.5*rhoMix*(u*u + v*v + w*w);
		consVec[iel*NEQ + 9] = 0; //vapor mass
		
	
		iel += gridDim.x * blockDim.x;
	}
}

void initialize(int nElem, double *cellCenter, double *consVec){
	
	initVec<<<(nElem + 511)/512,512>>>(d_input,nElem,input.NEQ,cellCenter,consVec);
		
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("CUDA error in initialize: %s\n", cudaGetErrorString(err));
	}
	
	err = cudaDeviceSynchronize();
	if (err != cudaSuccess) {
		    printf("KERNEL RUNTIME ERROR in initVec: %s\n", cudaGetErrorString(err));
	}

}
