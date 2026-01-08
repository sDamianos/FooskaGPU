
__device__ __forceinline__ void prim2cons(INPUT *d_input, float (&prim)[5], float *cons){

	cons[0] = prim[0];
	cons[1] = prim[0]*prim[1];
	cons[2] = prim[0]*prim[2];
	cons[3] = prim[0]*prim[3];
	cons[4] = prim[0]*(eqos(d_input,3,prim[4],prim[0]) + 0.5*( prim[1]*prim[1] + prim[2]*prim[2] + prim[3]*prim[3])); 

}
