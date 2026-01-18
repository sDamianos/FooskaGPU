
__device__ __forceinline__ void waveSpeeds(INPUT *d_input, float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&solL)[10], float (&solR)[10], float (&c)[3], float (&uL_lcl)[3], float (&uR_lcl)[3]){

	float cL,cR;
	float rhoL, rhoR;
	float pressL, pressR;
	float velL[3], velR[3];

	velL[0] = solL[3]; velL[1] = solL[4]; velL[2] = solL[5];
	velR[0] = solR[3]; velR[1] = solR[4]; velR[2] = solR[5];
	rotateVector_glb2lcl(velL,uL_lcl,n,nt1,nt2);
	rotateVector_glb2lcl(velR,uR_lcl,n,nt1,nt2);

	if((solL[8]+solR[8])/2 > d_input->chemLimit && d_input->woodsFlag == 1){
		float rhoMix, cLiq, cGas;
		rhoMix = solL[0]*solL[1] + (1-solL[0])*solL[2];
		cLiq = eqos(d_input,0,4,solL[1],solL[6]);
		cGas = eqos(d_input,1,4,solL[2],solL[7]);
		cL  = 1/sqrt(rhoMix*(solL[0]/(solL[1]*cLiq*cLiq) + (1-solL[0])/(solL[2]*cGas*cGas)));

		rhoMix = solR[0]*solR[1] + (1-solR[0])*solR[2];
		cLiq = eqos(d_input,0,4,solR[1],solR[6]);
		cGas = eqos(d_input,1,4,solR[2],solR[7]);
		cR  = 1/sqrt(rhoMix*(solR[0]/(solR[1]*cLiq*cLiq) + (1-solR[0])/(solR[2]*cGas*cGas)));
	}else{
		float YL,YG;
		YL = solL[0]*solL[1]/(solL[0]*solL[1] + (1-solL[0])*solL[2]);
		YG = (1-solL[0])*solL[2]/(solL[0]*solL[1] + (1-solL[0])*solL[2]);
		cL = sqrt(YL*pow(eqos(d_input,0,4,solL[1],solL[6]),2) + YG*pow(eqos(d_input,1,4,solL[2],solL[7]),2));
	
		YL = solR[0]*solR[1]/(solR[0]*solR[1] + (1-solR[0])*solR[2]);
		YG = (1-solR[0])*solR[2]/(solR[0]*solR[1] + (1-solR[0])*solR[2]);
		cR = sqrt(YL*pow(eqos(d_input,0,4,solR[1],solR[6]),2) + YG*pow(eqos(d_input,1,4,solR[2],solR[7]),2));
	}

	c[0] = fmin( uL_lcl[0] - cL , uR_lcl[0] - cR);   // Left  wave speed 
	c[2] = fmax( uL_lcl[0] + cL , uR_lcl[0] + cR);   // Right wave speed
    
    	pressL = solL[0]*solL[6] + (1-solL[0])*solL[7];
    	pressR = solR[0]*solR[6] + (1-solR[0])*solR[7]; 

	rhoL    = solL[0]*solL[1] + (1-solL[0])*solL[2];
	rhoR    = solR[0]*solR[1] + (1-solR[0])*solR[2]; 

	c[1] = ( pressR - pressL + rhoL*uL_lcl[0]*(c[0] - uL_lcl[0]) - rhoR*uR_lcl[0]*(c[2] - uR_lcl[0]) )/(rhoL*(c[0] - uL_lcl[0]) -  rhoR*(c[2] - uR_lcl[0]));	
}

