
__device__ __forceinline__ void bcInlet(INPUT *d_input, float (&solR)[5]){ //Dirichlet conditions

	float PressInlet = 4.2*101300; 
	float TempInlet  = 288;
//	float uInlet = 68.3;
//	float vInlet = 0.0;
//	float wInlet = 0.0;
	
	solR[0]  = eqos(d_input,2,PressInlet,TempInlet);
	float eInt = eqos(d_input,3,PressInlet,solR[0]);

//	solR[1] = uInlet;
//	solR[2] = vInlet;
//	solR[3] = wInlet;
	
	solR[4] = PressInlet;

}

__device__ __forceinline__ void bcOutlet(float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&solR)[5]){ // Neuman connditions (zero Gradient) so vecL == vecR and nothing has to change
	return;
}

__device__ __forceinline__ void bcWall(float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&solR)[5]){

	float uL_lcl[3],uR_lcl[3],solVel[3];
	
	solVel[0] = solR[1]; solVel[1] = solR[2]; solVel[2] = solR[3];
	rotateVector_glb2lcl(solVel,uL_lcl,n,nt1,nt2);

	uR_lcl[0] = -uL_lcl[0];
	uR_lcl[1] = -uL_lcl[1];
	uR_lcl[2] = -uL_lcl[2];
	
	rotateVector_lcl2glb(uR_lcl,solVel,n,nt1,nt2);
	solR[1] = solVel[0]; solR[2] = solVel[1]; solR[3] = solVel[2];

}

__device__ __forceinline__ void bcSlip(float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&solR)[5]){

	float uL_lcl[3],uR_lcl[3],solVel[3];
	
	solVel[0] = solR[1]; solVel[1] = solR[2]; solVel[2] = solR[3];
	rotateVector_glb2lcl(solVel,uL_lcl,n,nt1,nt2);

	uR_lcl[0] = -uL_lcl[0];
	uR_lcl[1] =  uL_lcl[1];
	uR_lcl[2] =  uL_lcl[2];

	rotateVector_lcl2glb(uR_lcl,solVel,n,nt1,nt2);
	solR[1] = solVel[0]; solR[2] = solVel[1]; solR[3] = solVel[2];
}

__device__ __forceinline__ void bcSymmetry(float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&solR)[5]){
	
	float uL_lcl[3],uR_lcl[3],solVel[3];
	
	solVel[0] = solR[1]; solVel[1] = solR[2]; solVel[2] = solR[3];
	rotateVector_glb2lcl(solVel,uL_lcl,n,nt1,nt2);

	uR_lcl[0] = -uL_lcl[0];
	uR_lcl[1] =  uL_lcl[1];
	uR_lcl[2] =  uL_lcl[2];

	rotateVector_lcl2glb(uR_lcl,solVel,n,nt1,nt2);
	solR[1] = solVel[0]; solR[2] = solVel[1]; solR[3] = solVel[2];


}
__device__ __forceinline__ void bcValues(INPUT *d_input,float (&n)[3], float (&nt1)[3], float (&nt2)[3], float (&solR)[5], int bc){
	if(bc==1){bcSymmetry(n,nt1,nt2,solR);}
	if(bc==2){bcOutlet(n,nt1,nt2,solR);}
	if(bc==3){bcInlet(d_input,solR);}
	if(bc==4){bcSlip(n,nt1,nt2,solR);}
	if(bc==5){bcWall(n,nt1,nt2,solR);}
}
