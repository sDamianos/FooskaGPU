// maybe i could predifine in the computeRHS finction the values rhoMix, pressMInx and carry them
// for flux computation and for star region. HOwever i am not sure how this will effect registers!
__device__ __forceinline__ void starRegion(INPUT *d_input, float (&sol)[10], float (&uLcl)[3], float SK,float S_S, float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&uNonCons)[3], float (&vecS)[10], float *consF)
{                   
	float u[3];
	float rhoSL,pSL,rhoSG,pSG;
	float press,rho,rhoS;

	u[0] = S_S; u[1] = uLcl[1]; u[2] = uLcl[2];
	rotateVector_lcl2glb(u,uNonCons,n,nt1,nt2);

	rhoSL    = sol[1] * (SK-uLcl[0])/(SK-S_S);
	rhoSG    = sol[2] * (SK-uLcl[0])/(SK-S_S);

	vecS[0] = sol[0];
	vecS[1] = sol[0]     * rhoSL;
	vecS[2] = (1-sol[0]) * rhoSG;

	rhoS = vecS[1] + vecS[2];
	vecS[3] = rhoS * uNonCons[0];
	vecS[4] = rhoS * uNonCons[1];
	vecS[5] = rhoS * uNonCons[2];
	
	pSL      = (sol[6] + d_input->pInf[0])*( (d_input->gamma[0]-1)*sol[1] - (d_input->gamma[0]+1)*rhoSL )/( (d_input->gamma[0]-1)*rhoSL - (d_input->gamma[0]+1)*sol[1] ) - d_input->pInf[0];	
	pSG      = (sol[7] + d_input->pInf[1])*( (d_input->gamma[1]-1)*sol[2] - (d_input->gamma[1]+1)*rhoSG )/( (d_input->gamma[1]-1)*rhoSG - (d_input->gamma[1]+1)*sol[2] ) - d_input->pInf[1];
	
	consF[6] = vecS[1] * eqos(d_input,0,3,pSL,rhoSL);
	consF[7] = vecS[2] * eqos(d_input,1,3,pSG,rhoSG);

	press   = sol[0] * sol[6] + (1-sol[0])* sol[7];
	rho     = sol[0] * sol[1] + (1-sol[0])* sol[2];
	vecS[8] = ( (SK-uLcl[0])/(SK-S_S) ) * ( consF[8] + (S_S-uLcl[0])*( S_S + press/(rho*(SK-uLcl[0])) ) );
	
	vecS[9] = sol[8] * (SK-uLcl[0])/(SK-S_S); 
}
