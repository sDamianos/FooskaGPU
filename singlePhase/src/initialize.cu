#include "strdata.h"
#include "eqos.h"

__constant__ float press1 = 4.2*101300;     __constant__ float press2 = 101300;
__constant__ float temp1  = 288;        __constant__ float temp2  = 288;
__constant__ float u1     = 0.0;        __constant__ float u2    = 0.0;
__constant__ float v1     = 0.0;        __constant__ float v2     = 0.0;
__constant__ float w1     = 0.0;        __constant__ float w2     = 0.0;

__global__ void initVec(INPUT *d_input, int nElem, float *cellCenter, float *consVec){
	
	int iel = threadIdx.x + blockIdx.x * blockDim.x;
	float press,temp,rho,eInt,u,v,w;
	int is_region_1;
	
	while(iel < nElem){
		
		is_region_1 = (cellCenter[iel*3 + 0] < 0.0); 
		
		press = is_region_1 ? press1 : press2;
		temp  = is_region_1 ? temp1  : temp2;
		u     = is_region_1 ? u1     : u2;
		v     = is_region_1 ? v1     : v2;
		w     = is_region_1 ? w1     : w2;

		rho  = eqos(d_input,2,press,temp);
		eInt = eqos(d_input,3,press,rho);

		consVec[iel*5 + 0] = rho;
		consVec[iel*5 + 1] = rho*u;
		consVec[iel*5 + 2] = rho*v;
		consVec[iel*5 + 3] = rho*w;
		consVec[iel*5 + 4] = rho*(eInt + 0.5*u*u + v*v + w*w);
		
		iel += gridDim.x * blockDim.x;
	}
}
	

// CPU version its actual to fast, (3 orders of magnitude lower), and also this step
// is fully scalable. However also teh GPU step is fast (~1-2 ms) so we could say that
// is irrelevant what we wil use.
void initialize(int nElem, float *cellCenter, float *consVec){

	initVec<<<(nElem + 511)/512,512>>>(d_input,nElem,cellCenter,consVec);

}
/*	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);	
	cudaEventRecord(start,0);	
	
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float   elapsedTime;
	cudaEventElapsedTime(&elapsedTime,start,stop);
	printf( "Time to generate:  %3.5f ms\n", elapsedTime );
*/
