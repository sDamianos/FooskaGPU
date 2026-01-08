
__device__ __forceinline__ void cons2prim(INPUT *d_input, float *consVec, float (&sol)[5]){
			
	sol[0] = consVec[0];
	sol[1] = consVec[1]/sol[0];
	sol[2] = consVec[2]/sol[0];
	sol[3] = consVec[3]/sol[0];
	sol[4] = consVec[4]/sol[0] - 0.5*(sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3]);
	sol[4] = eqos(d_input,0, sol[0],sol[4]);
}
