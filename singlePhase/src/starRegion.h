__device__ __forceinline__ void starRegion(INPUT *d_input, float (&sol)[5], float (&uLcl)[3], float SK, float S_S, float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&vecS)[5]){
                   
	float u[3],vecGlb[3];

	u[0] = S_S; u[1] = uLcl[1]; u[2] = uLcl[2];
	rotateVector_lcl2glb(u,vecGlb,n,nt1,nt2);

	
	vecS[0] = sol[0]  * (SK-uLcl[0])/(SK-S_S);           
	vecS[1] = vecS[0] * vecGlb[0];
	vecS[2] = vecS[0] * vecGlb[1];
	vecS[3] = vecS[0] * vecGlb[2];

	vecS[4] = vecS[0]* ( eqos(d_input,3,sol[4],sol[0]) + 0.5*(sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3]) + (S_S-uLcl[0])*(S_S + sol[4]/(sol[0]*(SK-uLcl[0]))) );	
}
