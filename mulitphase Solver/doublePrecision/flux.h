__device__ __forceinline__ void fluxSource(INPUT *d_input, double *uNonCons, double *consVecL, double *consVecR, double (&n)[3], double (&fluxNonCons)[6]){

	double avf,rhoLiq,rhoGas,eLiq,eGas,pLiq,pGas;

	// Left State
	avf    = fmin(fmax(d_input->avfMin/100,consVecL[0]),1-d_input->avfMin/100);
	rhoLiq = consVecL[1]/avf;
	rhoGas = consVecL[2]/(1-avf);
	eLiq   = consVecL[6]/(avf*rhoLiq);
	eGas   = consVecL[7]/((1-avf)*rhoGas);
	pLiq   = eqos(d_input,0,0,rhoLiq,eLiq);
	pGas   = eqos(d_input,1,0,rhoGas,eGas);

	fluxNonCons[0*2 + 0] = avf*uNonCons[0]*n[0]          + avf*uNonCons[1]*n[1]          + avf*uNonCons[2]*n[2];
	fluxNonCons[1*2 + 0] = avf*pLiq*uNonCons[0]*n[0]     + avf*pLiq*uNonCons[1]*n[1]     + avf*pLiq*uNonCons[2]*n[2];
	fluxNonCons[2*2 + 0] = (1-avf)*pGas*uNonCons[0]*n[0] + (1-avf)*pGas*uNonCons[1]*n[1] + (1-avf)*pGas*uNonCons[2]*n[2];
	
	// Right State
	avf    = fmin(fmax(d_input->avfMin/100,consVecR[0]),1-d_input->avfMin/100);
	rhoLiq = consVecR[1]/avf;
	rhoGas = consVecR[2]/(1-avf);
	eLiq   = consVecR[6]/(avf*rhoLiq);
	eGas   = consVecR[7]/((1-avf)*rhoGas);
	pLiq   = eqos(d_input,0,0,rhoLiq,eLiq);
	pGas   = eqos(d_input,1,0,rhoGas,eGas);
	
	fluxNonCons[0*2 + 1] = avf*uNonCons[0]*n[0]          + avf*uNonCons[1]*n[1]          + avf*uNonCons[2]*n[2];
	fluxNonCons[1*2 + 1] = avf*pLiq*uNonCons[0]*n[0]     + avf*pLiq*uNonCons[1]*n[1]     + avf*pLiq*uNonCons[2]*n[2];
	fluxNonCons[2*2 + 1] = (1-avf)*pGas*uNonCons[0]*n[0] + (1-avf)*pGas*uNonCons[1]*n[1] + (1-avf)*pGas*uNonCons[2]*n[2];
}

__device__ __forceinline__ void flux(INPUT *d_input, int isStarRegion, double (&primVec)[10], double (&fluxAdv)[15], double (&n)[3], double (&uNonCons)[3], double areL, double areG){

	double rho, press, rE;
		
	rho   = primVec[0]*primVec[1] + (1-primVec[0])*primVec[2];
	press = primVec[0]*primVec[6] + (1-primVec[0])*primVec[7];
	rE = primVec[0]*primVec[1]*eqos(d_input,0,3,primVec[6],primVec[1]) + 
	    (1-primVec[0])*primVec[2]*eqos(d_input,1,3,primVec[7],primVec[2]) + 
	    0.5*rho*( primVec[3]*primVec[3] + primVec[4]*primVec[4] + primVec[5]*primVec[5] );

//	volume fraction  
	fluxAdv[0]  = uNonCons[0]*primVec[0]*n[0];
	fluxAdv[0] += uNonCons[1]*primVec[0]*n[1];
	fluxAdv[0] += uNonCons[2]*primVec[0]*n[2];
	
//	Mass Liquid
	fluxAdv[1]  = primVec[3]*primVec[0]*primVec[1]*n[0];
	fluxAdv[1] += primVec[4]*primVec[0]*primVec[1]*n[1];
	fluxAdv[1] += primVec[5]*primVec[0]*primVec[1]*n[2];

//	Mass Gas
	fluxAdv[2]  = primVec[3]*(1-primVec[0])*primVec[2]*n[0];
	fluxAdv[2] += primVec[4]*(1-primVec[0])*primVec[2]*n[1];
	fluxAdv[2] += primVec[5]*(1-primVec[0])*primVec[2]*n[2];

//	Momentum X
	fluxAdv[3]  = rho*primVec[3]*primVec[3]*n[0];
	fluxAdv[3] += rho*primVec[3]*primVec[4]*n[1];
	fluxAdv[3] += rho*primVec[3]*primVec[5]*n[2];

//	Momentum Y
	fluxAdv[4]  = rho*primVec[4]*primVec[3]*n[0];
	fluxAdv[4] += rho*primVec[4]*primVec[4]*n[1];
	fluxAdv[4] += rho*primVec[4]*primVec[5]*n[2];

//	Momentum Z
	fluxAdv[5]  = rho*primVec[5]*primVec[3]*n[0];
	fluxAdv[5] += rho*primVec[5]*primVec[4]*n[1];
	fluxAdv[5] += rho*primVec[5]*primVec[5]*n[2];

//	Liquid Phasic Energy
	fluxAdv[6]  = uNonCons[0]*areL*n[0];
	fluxAdv[6] += uNonCons[1]*areL*n[1];
	fluxAdv[6] += uNonCons[2]*areL*n[2];

//	Gas Phasic Energy
	fluxAdv[7]  = uNonCons[0]*areG*n[0];
	fluxAdv[7] += uNonCons[1]*areG*n[1];
	fluxAdv[7] += uNonCons[2]*areG*n[2];

//	Total Energy
	fluxAdv[8]  = primVec[3]*(rE + press)*n[0];
	fluxAdv[8] += primVec[4]*(rE + press)*n[1];
	fluxAdv[8] += primVec[5]*(rE + press)*n[2];

//	vapor Mass
	fluxAdv[9]  = primVec[3]*primVec[8]*n[0];
	fluxAdv[9] += primVec[4]*primVec[8]*n[1];
	fluxAdv[9] += primVec[5]*primVec[8]*n[2];	
}

/*
//	volume fraction  
	fluxAdv[0*3 + 0] = uNonCons[0]*(primVec[0]-sTerm[0])*n[0];
	fluxAdv[0*3 + 1] = uNonCons[1]*(primVec[0]-sTerm[0])*n[1];
	fluxAdv[0*3 + 2] = uNonCons[2]*(primVec[0]-sTerm[0])*n[2];
	
//	Mass Liquid
	fluxAdv[1*3 + 0] = primVec[3]*primVec[0]*primVec[1]*n[0];
	fluxAdv[1*3 + 1] = primVec[4]*primVec[0]*primVec[1]*n[1];
	fluxAdv[1*3 + 2] = primVec[5]*primVec[0]*primVec[1]*n[2];

//	Mass Gas
	fluxAdv[2*3 + 0] = primVec[3]*(1-primVec[0])*primVec[2]*n[0];
	fluxAdv[2*3 + 1] = primVec[4]*(1-primVec[0])*primVec[2]*n[1];
	fluxAdv[2*3 + 2] = primVec[5]*(1-primVec[0])*primVec[2]*n[2];

//	Momentum X
	fluxAdv[3*3 + 0] = (rho*primVec[3]*primVec[3] + press)*n[0];
	fluxAdv[3*3 + 1] =  rho*primVec[3]*primVec[4]*n[1];
	fluxAdv[3*3 + 2] =  rho*primVec[3]*primVec[5]*n[2];

//	Momentum Y
	fluxAdv[4*3 + 0] =  rho*primVec[2]*primVec[1]*n[0];
	fluxAdv[4*3 + 1] = (rho*primVec[2]*primVec[2] + press)*n[1];
	fluxAdv[4*3 + 2] =  rho*primVec[2]*primVec[3]*n[2];

//	Momentum Z
	fluxAdv[5*3 + 0] =  rho*primVec[3]*primVec[1]*n[0];
	fluxAdv[5*3 + 1] =  rho*primVec[3]*primVec[2]*n[1];
	fluxAdv[5*3 + 2] = (rho*primVec[3]*primVec[3] + press)*n[2];

//	Liquid Phasic Energy
	fluxAdv[6*3 + 0] = uNonCons[0]*(..)*n[0];
	fluxAdv[6*3 + 1] = uNonCons[1]*(..)*n1;
	fluxAdv[6*3 + 2] = uNonCons[2]*(..)*n2;

//	Liquid Phasic Energy
	fluxAdv[7*3 + 0] = uNonCons[0]*(..)*n0;
	fluxAdv[7*3 + 1] = uNonCons[1]*(..)*n1;
	fluxAdv[7*3 + 2] = uNonCons[2]*(..)*n2;

//	Total Energy
	fluxAdv[8*3 + 0] = n[0]*primVec[3]*(rE + press);
	fluxAdv[8*3 + 1] = n[1]*primVec[4]*(rE + press);
	fluxAdv[8*3 + 2] = n[2]*primVec[5]*(rE + press);

//	vapor Mass
	fluxAdv[9*3 + 0] = primVec[3]*primVec[8]*n[0];
	fluxAdv[9*3 + 1] = primVec[4]*primVec[8]*n[1];
	fluxAdv[9*3 + 2] = primVec[5]*primVec[8]*n[2];
*/
