__device__ __forceinline__ void bcInlet(INPUT *d_input, double (&solR)[10]){ //Dirichlet conditions

	
 	float avfInlet   = 1.0 - d_input->avfMin;
	float PressInlet = 11.4e5; 
	float TempInlet  = 115.3;
//	float uInlet = 68.3;
//	float vInlet = 0.0;
//	float wInlet = 0.0;
	
	solR[0]  = avfInlet;
	solR[1]  = eqos(d_input,0,2,PressInlet,TempInlet);
	solR[2]  = eqos(d_input,1,2,PressInlet,TempInlet);

//	solR[3] = uInlet;
//	solR[4] = vInlet;
//	solR[5] = wInlet;
	
	solR[6] = PressInlet;
	solR[7] = PressInlet;
	solR[8] = 0;
	solR[9] = 0;

}

__device__ __forceinline__ void bcOutlet(float (&n)[3], float (&nt1)[3], float (&nt2)[3], double (&solR)[10]){ // Neuman connditions (zero Gradient) so vecL == vecR and nothing has to change
	return;
}

__device__ __forceinline__ void bcWall(float (&n)[3], float (&nt1)[3], float (&nt2)[3], double (&solR)[10]){

	float uL_lcl[3],uR_lcl[3],solVel[3];
	
	solVel[0] = solR[3]; solVel[1] = solR[4]; solVel[2] = solR[5];
	rotateVector_glb2lcl(solVel,uL_lcl,n,nt1,nt2);

	uR_lcl[0] = -uL_lcl[0];
	uR_lcl[1] = -uL_lcl[1];
	uR_lcl[2] = -uL_lcl[2];
	
	rotateVector_lcl2glb(uR_lcl,solVel,n,nt1,nt2);
	solR[3] = solVel[0]; solR[4] = solVel[1]; solR[5] = solVel[2];

}

__device__ __forceinline__ void bcSlip(float (&n)[3], float (&nt1)[3], float (&nt2)[3], double (&solR)[10]){

	float uL_lcl[3],uR_lcl[3],solVel[3];
	
	solVel[0] = solR[3]; solVel[1] = solR[4]; solVel[2] = solR[5];
	rotateVector_glb2lcl(solVel,uL_lcl,n,nt1,nt2);

	uR_lcl[0] = -uL_lcl[0];
	uR_lcl[1] =  uL_lcl[1];
	uR_lcl[2] =  uL_lcl[2];

	rotateVector_lcl2glb(uR_lcl,solVel,n,nt1,nt2);
	solR[3] = solVel[0]; solR[4]= solVel[1]; solR[5] = solVel[2];
}

__device__ __forceinline__ void bcSymmetry(float (&n)[3], float (&nt1)[3], float (&nt2)[3], double (&solR)[10]){
	
	float uL_lcl[3],uR_lcl[3],solVel[3];
	
	solVel[0] = solR[3]; solVel[1] = solR[4]; solVel[2] = solR[5];
	rotateVector_glb2lcl(solVel,uL_lcl,n,nt1,nt2);

	uR_lcl[0] = -uL_lcl[0];
	uR_lcl[1] =  uL_lcl[1];
	uR_lcl[2] =  uL_lcl[2];

	rotateVector_lcl2glb(uR_lcl,solVel,n,nt1,nt2);
	solR[3] = solVel[0]; solR[4] = solVel[1]; solR[5] = solVel[2];


}

__device__ __forceinline__ void bcValues(INPUT *d_input,float (&n)[3], float (&nt1)[3], float (&nt2)[3], double (&solR)[10], int bc){
	if(bc==1){bcSymmetry(n,nt1,nt2,solR);}
	if(bc==2){bcOutlet(n,nt1,nt2,solR);}
	if(bc==3){bcInlet(d_input,solR);}
	if(bc==4){bcSlip(n,nt1,nt2,solR);}
	if(bc==5){bcWall(n,nt1,nt2,solR);}
}

