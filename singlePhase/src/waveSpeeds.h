__device__ __forceinline__ void waveSpeeds(INPUT *d_input, float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&solL)[5], float (&solR)[5], float (&c)[3], float (&uL_lcl)[3], float (&uR_lcl)[3]){

	float cL,cR;
	float velL[3], velR[3];

	velL[0] = solL[1]; velL[1] = solL[2]; velL[2] = solL[3];
	velR[0] = solR[1]; velR[1] = solR[2]; velR[2] = solR[3];
	rotateVector_glb2lcl(velL,uL_lcl,n,nt1,nt2);
	rotateVector_glb2lcl(velR,uR_lcl,n,nt1,nt2);

	cL = eqos(d_input,4,solL[0],solL[4]);
	cR = eqos(d_input,4,solR[0],solR[4]);

	c[0] = fmin( uL_lcl[0] - cL , uR_lcl[0] - cR);   // Left  wave speed 
	c[2] = fmax( uL_lcl[0] + cL , uR_lcl[0] + cR);   // Right wave speed
    
	c[1] = ( solR[4] - solL[4] + solL[0]*uL_lcl[0]*(c[0] - uL_lcl[0]) - solR[0]*uR_lcl[0]*(c[2] - uR_lcl[0]) )/(solL[0]*(c[0] - uL_lcl[0]) -  solR[0]*(c[2] - uR_lcl[0]));	
}

