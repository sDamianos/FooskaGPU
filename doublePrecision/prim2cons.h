
__device__ __forceinline__ void prim2cons(INPUT *d_input, double (&prim)[10], double *cons){

	double rho;

	cons[0] = prim[0];
	cons[1] = prim[0]*prim[1];
	cons[2] = (1-prim[0])*prim[2];

	rho = cons[1] + cons[2];
	cons[3] = rho*prim[3];
	cons[4] = rho*prim[4];
	cons[5] = rho*prim[5];

	cons[6] = cons[1]*eqos(d_input,0,3,prim[6],prim[1]);
	cons[7] = cons[2]*eqos(d_input,1,3,prim[7],prim[2]);

	cons[8] = cons[6] + cons[7] + 0.5*rho*(prim[3]*prim[3] + prim[4]*prim[4] + prim[5]*prim[5]);
	cons[9] = prim[8];

}

