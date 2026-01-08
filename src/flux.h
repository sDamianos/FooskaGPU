
__device__ __forceinline__ void flux(INPUT *d_input,float (&primVec)[5], float (&fluxAdv)[15], float (&n)[3]) {
  
	fluxAdv[0*3 + 0] = primVec[1]*primVec[0]*n[0];
	fluxAdv[0*3 + 1] = primVec[2]*primVec[0]*n[1];
	fluxAdv[0*3 + 2] = primVec[3]*primVec[0]*n[2];
	
	fluxAdv[1*3 + 0] = (primVec[0]*primVec[1]*primVec[1] + primVec[4])*n[0];
	fluxAdv[1*3 + 1] = primVec[0]*primVec[1]*primVec[2]*n[1];
	fluxAdv[1*3 + 2] = primVec[0]*primVec[1]*primVec[3]*n[2];

	fluxAdv[2*3 + 0] = primVec[0]*primVec[2]*primVec[1]*n[0];
	fluxAdv[2*3 + 1] = (primVec[0]*primVec[2]*primVec[2] + primVec[4])*n[1];
	fluxAdv[2*3 + 2] = primVec[0]*primVec[2]*primVec[3]*n[2];

	fluxAdv[3*3 + 0] = primVec[0]*primVec[3]*primVec[1]*n[0];
	fluxAdv[3*3 + 1] = primVec[0]*primVec[3]*primVec[2]*n[1];
	fluxAdv[3*3 + 2] = (primVec[0]*primVec[3]*primVec[3] + primVec[4])*n[2];

	fluxAdv[4*3 + 0] = n[0]*primVec[1]*(primVec[0]*( eqos(d_input,3,primVec[4],primVec[0]) + 0.5*( primVec[1]*primVec[1] + primVec[2]*primVec[2] + primVec[3]*primVec[3] ) ) + primVec[4]);
	fluxAdv[4*3 + 1] = n[1]*primVec[2]*(primVec[0]*( eqos(d_input,3,primVec[4],primVec[0]) + 0.5*( primVec[1]*primVec[1] + primVec[2]*primVec[2] + primVec[3]*primVec[3] ) ) + primVec[4]);
	fluxAdv[4*3 + 2] = n[2]*primVec[3]*(primVec[0]*( eqos(d_input,3,primVec[4],primVec[0]) + 0.5*( primVec[1]*primVec[1] + primVec[2]*primVec[2] + primVec[3]*primVec[3] ) ) + primVec[4]);
}
