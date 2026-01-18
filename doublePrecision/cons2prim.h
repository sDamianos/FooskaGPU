
__device__ __forceinline__ void cons2prim(INPUT *d_input, double *consVec, double (&sol)[10]){
			

	double rhoMix, eIntL, eIntG;

	sol[0] = fmin(fmax(d_input->avfMin/100,consVec[0]),1.0-d_input->avfMin/100);
	sol[1] = consVec[1]/sol[0];
	sol[2] = consVec[2]/(1-sol[0]);
	
	rhoMix = sol[0]*sol[1] + (1-sol[0])*sol[2];
	sol[3] = consVec[3]/rhoMix;
	sol[4] = consVec[4]/rhoMix;
	sol[5] = consVec[5]/rhoMix;

	eIntL  = consVec[6]/(sol[0]*sol[1]);
	eIntG  = consVec[7]/((1-sol[0])*sol[2]); 
	sol[6] = eqos(d_input,0,0,sol[1],eIntL);
	sol[7] = eqos(d_input,1,0,sol[2],eIntG);

	sol[8] = consVec[9];
	sol[9] = rhoMix; //or Temp could also work
//	sol[9] = sol[0] * eqos(d_input,0,1,sol[1],eIntL) + sol[1] * eqos(d_input,1,1,sol[1],eIntL)
}

