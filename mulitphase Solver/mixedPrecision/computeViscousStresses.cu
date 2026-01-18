#include "strdata.h"

__device__ __forceinline__  float computeEddyViscocity(float ux, float uy, float uz, float vx, float vy, float vz, float wx, float wy, float wz, float volume){

	double g[3][3], square_g[3][3], s_d[3][3], strain[3][3];
	double T, num, denom, Ls, sum1, sum2;
	double eddyVisc;
	int kron;

	g[0][0] = (double) ux; g[0][1] = (double) uy; g[0][2] = (double) uz;
	g[1][0] = (double) vx; g[1][1] = (double) vy; g[1][2] = (double) vz;
	g[2][0] = (double) wx; g[2][1] = (double) wy; g[2][2] = (double) wz;
        
        for(int i = 0; i < 3; i++) {
                for(int j = 0; j < 3; j++) {
                        square_g[i][j] = 0.0;
                        for(int k = 0; k < 3; k++) {
                                square_g[i][j] += g[i][k] * g[k][j];
                        }
                }
        }
	
	T = 0.0;
	for(int i = 0; i < 3; i++){
		T += square_g[i][i];
	}
 
 	kron = 0;
	for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; ++j) {
			kron = ( i == j);
                        s_d[i][j] = 0.5*(square_g[i][j] + square_g[j][i])-(1.0/3.0)*kron*T;
                }
        }
	
        for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; j++){
                        strain[i][j] = 0.5*(g[i][j] + g[j][i]);
                }
        }
	
	sum1 = 0;
	sum2 = 0; 
        for(int i = 0; i < 3; i++){
                for(int j = 0; j < 3; ++j) {
                        sum1 += s_d[i][j]*s_d[i][j];
                        sum2 += strain[i][j]*strain[i][j];
                }
        }
	num = pow(sum1,1.5);
	denom = pow(sum1,1.25) + pow (sum2,2.5) + 1e-12;
	if(denom < 1e-10){ num = 0;}

        Ls = 0.325*pow((double)volume,0.3333333333333);
        
	eddyVisc = Ls*Ls*num/denom;

	return eddyVisc;

}

// Prallelizing it per iDirections increase the computational cost by 25% and doesn't decrease the number of registers.
// Basically i have still to define the whole stress tensor vector, so the only thing that i am saving is 
// the loop over the lines of teh viscocity matix. This is just a small fraction of the computational 
// cost and thus is not really scalable. However I could also try to spot if there are differencies when runnning to A100.
__global__ void computeViscousStresses(INPUT *d_input, int nFcs, int *fc2el, int *boundCond, double *primVecF, float *grad, float *n, float *area, float *volume, float *RHS){
	
	int ifc = threadIdx.x + blockIdx.x * blockDim.x;

	float stress_tens[3][3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	float visc_stress[3][3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	float nVec[3];
	float fluxViscMom[3];
	float velF[3];
	float stress_times_vel[3];
	float visc_stress_en_tot[3];
	float fluxViscEner[3];
	float ux,uy,uz,vx,vy,vz,wx,wy,wz;
	float avf1,avf2,rhoL,rhoG,rho1,rho2;
	float rhoF,avfF;
	float eddyVisc=0;
	float visc;
	float areaF;
	int id =0;
	int bc;

	while(ifc < nFcs){
		bc = boundCond[ifc];
		id = (bc == 0);
		
		ux = (grad[fc2el[ifc*2+0]*10*3 + 3*3 + 0] + grad[fc2el[ifc*2+id]*10*3 + 3*3 + 0])/2;
		uy = (grad[fc2el[ifc*2+0]*10*3 + 3*3 + 1] + grad[fc2el[ifc*2+id]*10*3 + 3*3 + 1])/2;
		uz = (grad[fc2el[ifc*2+0]*10*3 + 3*3 + 2] + grad[fc2el[ifc*2+id]*10*3 + 3*3 + 2])/2;

		vx = (grad[fc2el[ifc*2+0]*10*3 + 4*3 + 0] + grad[fc2el[ifc*2+id]*10*3 + 4*3 + 0])/2;
		vy = (grad[fc2el[ifc*2+0]*10*3 + 4*3 + 1] + grad[fc2el[ifc*2+id]*10*3 + 4*3 + 1])/2;
		vz = (grad[fc2el[ifc*2+0]*10*3 + 4*3 + 2] + grad[fc2el[ifc*2+id]*10*3 + 4*3 + 2])/2;
		
		wx = (grad[fc2el[ifc*2+0]*10*3 + 5*3 + 0] + grad[fc2el[ifc*2+id]*10*3 + 5*3 + 0])/2;
		wy = (grad[fc2el[ifc*2+0]*10*3 + 5*3 + 1] + grad[fc2el[ifc*2+id]*10*3 + 5*3 + 1])/2;
		wz = (grad[fc2el[ifc*2+0]*10*3 + 5*3 + 2] + grad[fc2el[ifc*2+id]*10*3 + 5*3 + 2])/2;	

		if(d_input->turbulence == 1){ eddyVisc = (float) computeEddyViscocity(ux,uy,uz,vx,vy,vz,wx,wy,wz,(volume[fc2el[ifc*2+0]]+volume[fc2el[ifc*2+id]])/2);}

		stress_tens[0][0] = ( 2.0*ux - (2.0/3.0)* (ux+vy+wz) );
		stress_tens[1][1] = ( 2.0*vy - (2.0/3.0)* (ux+vy+wz) );
		stress_tens[2][2] = ( 2.0*wz - (2.0/3.0)* (ux+vy+wz) );

		stress_tens[0][1] = (uy+vx);
		stress_tens[0][2] = (uz+wx);
		stress_tens[1][2] = (vz+wy);

		stress_tens[1][0] = stress_tens[0][1];
		stress_tens[2][1] = stress_tens[1][2];
		stress_tens[2][0] = stress_tens[0][2];
		
		avf1 = primVecF[ifc*2*10 + 0*10 + 0];
		rhoL = primVecF[ifc*2*10 + 0*10 + 1];
		rhoG = primVecF[ifc*2*10 + 0*10 + 2];
		rho1 = avf1* rhoL + (1-avf1) * rhoG;
		
		avf2 = primVecF[ifc*2*10 + 1*10 + 0];
		rhoL = primVecF[ifc*2*10 + 1*10 + 1];
		rhoG = primVecF[ifc*2*10 + 1*10 + 2];
		rho2 = avf2 * rhoL + (1-avf2) * rhoG;
	
		rhoF = (rho1+rho2)/2;
		avfF = (avf1+avf2)/2;
		visc = avfF * d_input->dynVisc[0] + (1-avfF) * d_input->dynVisc[1] + rhoF*eddyVisc;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				visc_stress[i][j] = visc * stress_tens[i][j];
			}
		}

		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];

		fluxViscMom[0] =  visc_stress[0][0]*nVec[0] + visc_stress[0][1]*nVec[1] + visc_stress[0][2]*nVec[2];
		fluxViscMom[1] =  visc_stress[1][0]*nVec[0] + visc_stress[1][1]*nVec[1] + visc_stress[1][2]*nVec[2];
		fluxViscMom[2] =  visc_stress[2][0]*nVec[0] + visc_stress[2][1]*nVec[1] + visc_stress[2][2]*nVec[2];
	
		areaF = area[ifc];	
		for(int iv = 3; iv < 6; iv++){	
			atomicAdd( &RHS[fc2el[ifc*2+0]*10+iv], (fluxViscMom[iv-3])*areaF);
			if(bc==0){atomicAdd( &RHS[fc2el[ifc*2+1]*10+iv], -(fluxViscMom[iv-3])*areaF);}
		}
		
		velF[0] = (primVecF[ifc*2*10 + 0*10 + 3] + primVecF[ifc*2*10 + 1*10 + 3])/2;
		velF[1] = (primVecF[ifc*2*10 + 0*10 + 4] + primVecF[ifc*2*10 + 1*10 + 4])/2;
		velF[2] = (primVecF[ifc*2*10 + 0*10 + 5] + primVecF[ifc*2*10 + 1*10 + 5])/2;
			
		for (int i = 0; i < 3; i++) {
			stress_times_vel[i] = 0;
			for (int j = 0; j < 3; j++) {
				stress_times_vel[i] += stress_tens[i][j] * velF[j];
			}
			visc_stress_en_tot[i] = visc * stress_times_vel[i];
		}

		fluxViscEner[0] = visc_stress_en_tot[0]*nVec[0];
		fluxViscEner[1] = visc_stress_en_tot[1]*nVec[1];
		fluxViscEner[2] = visc_stress_en_tot[2]*nVec[2];
		
		bc = boundCond[ifc];
		areaF = area[ifc];
		atomicAdd( &RHS[fc2el[ifc*2+0]*10+8], (fluxViscEner[0] + fluxViscEner[1] + fluxViscEner[2])*areaF);
		if(bc==0){atomicAdd( &RHS[fc2el[ifc*2+1]*10+8], -(fluxViscEner[0] + fluxViscEner[1] + fluxViscEner[2])*areaF);}

		ifc += blockDim.x * gridDim.x;
	}
}


__global__ void computeViscousStressesProx(INPUT *d_input, int start, int nFcs, int max, int *fc2el, int *neigRank4fc, int *lclFc2idRcv, double *primVecF, float *grad, float *n, float *area, float *volume, float *RHS, float *recvBuffGrad, double *recvBuffPrim){
	
	int ifc = start + threadIdx.x + blockIdx.x * blockDim.x;

	float stress_tens[3][3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	float visc_stress[3][3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	float nVec[3];
	float fluxViscMom[3];
	float velF[3];
	float stress_times_vel[3];
	float visc_stress_en_tot[3];
	float fluxViscEner[3];
	float ux,uy,uz,vx,vy,vz,wx,wy,wz;
	float avf1,avf2,rhoL,rhoG,rho1,rho2;
	float rhoF,avfF;
	float eddyVisc=0;
	float visc;
	float areaF;
	int neigRank;
	int iProx;

	while(ifc < nFcs){
		
		neigRank = neigRank4fc[ifc];
		iProx    = lclFc2idRcv[ifc];
		
		ux = (grad[fc2el[ifc*2+0]*10*3 + 3*3 + 0] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 0*3 + 0])/2;
		uy = (grad[fc2el[ifc*2+0]*10*3 + 3*3 + 1] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 0*3 + 1])/2;
		uz = (grad[fc2el[ifc*2+0]*10*3 + 3*3 + 2] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 0*3 + 2])/2;

		vx = (grad[fc2el[ifc*2+0]*10*3 + 4*3 + 0] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 1*3 + 0])/2;
		vy = (grad[fc2el[ifc*2+0]*10*3 + 4*3 + 1] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 1*3 + 1])/2;
		vz = (grad[fc2el[ifc*2+0]*10*3 + 4*3 + 2] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 1*3 + 2])/2;
	
		wx = (grad[fc2el[ifc*2+0]*10*3 + 5*3 + 0] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 2*3 + 0])/2;
		wy = (grad[fc2el[ifc*2+0]*10*3 + 5*3 + 1] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 2*3 + 1])/2;
		wz = (grad[fc2el[ifc*2+0]*10*3 + 5*3 + 2] + recvBuffGrad[neigRank*max*3*3 + iProx*3*3 + 2*3 + 2])/2;	

		if(d_input->turbulence == 1){ eddyVisc = (float) computeEddyViscocity(ux,uy,uz,vx,vy,vz,wx,wy,wz,(volume[fc2el[ifc*2+0]]+volume[fc2el[ifc*2+1]])/2);}

		stress_tens[0][0] = ( 2.0*ux - (2.0/3.0)* (ux+vy+wz) );
		stress_tens[1][1] = ( 2.0*vy - (2.0/3.0)* (ux+vy+wz) );
		stress_tens[2][2] = ( 2.0*wz - (2.0/3.0)* (ux+vy+wz) );

		stress_tens[0][1] = (uy+vx);
		stress_tens[0][2] = (uz+wx);
		stress_tens[1][2] = (vz+wy);

		stress_tens[1][0] = stress_tens[0][1];
		stress_tens[2][1] = stress_tens[1][2];
		stress_tens[2][0] = stress_tens[0][2];
	
		avf1 = primVecF[ifc*2*10 + 0*10 + 0];
		rhoL = primVecF[ifc*2*10 + 0*10 + 1];
		rhoG = primVecF[ifc*2*10 + 0*10 + 2];
		rho1 = avf1* rhoL + (1-avf1) * rhoG;
		
		avf2 = recvBuffPrim[neigRank*max*10 + iProx*10 + 0];
		rhoL = recvBuffPrim[neigRank*max*10 + iProx*10 + 1];
		rhoG = recvBuffPrim[neigRank*max*10 + iProx*10 + 2]; 
		rho2 = avf2 * rhoL + (1-avf2) * rhoG;
	
		rhoF = (rho1+rho2)/2;
		avfF = (avf1+avf2)/2;
		visc = avfF * d_input->dynVisc[0] + (1-avfF) * d_input->dynVisc[1] + rhoF*eddyVisc;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				visc_stress[i][j] = visc * stress_tens[i][j];
			}
		}

		nVec[0] = n[ifc*3 + 0];
		nVec[1] = n[ifc*3 + 1];
		nVec[2] = n[ifc*3 + 2];

		fluxViscMom[0] =  visc_stress[0][0]*nVec[0] + visc_stress[0][1]*nVec[1] + visc_stress[0][2]*nVec[2];
		fluxViscMom[1] =  visc_stress[1][0]*nVec[0] + visc_stress[1][1]*nVec[1] + visc_stress[1][2]*nVec[2];
		fluxViscMom[2] =  visc_stress[2][0]*nVec[0] + visc_stress[2][1]*nVec[1] + visc_stress[2][2]*nVec[2];
	
		areaF = area[ifc];	
		for(int iv = 3; iv < 6; iv++){	
			atomicAdd( &RHS[fc2el[ifc*2+0]*10+iv], (fluxViscMom[iv-3])*areaF);
		}
		
		velF[0] = (primVecF[ifc*2*10 + 0*10 + 3] + recvBuffPrim[neigRank*max*10 + iProx*10 + 3])/2;
		velF[1] = (primVecF[ifc*2*10 + 0*10 + 4] + recvBuffPrim[neigRank*max*10 + iProx*10 + 4])/2;
		velF[2] = (primVecF[ifc*2*10 + 0*10 + 5] + recvBuffPrim[neigRank*max*10 + iProx*10 + 5])/2;
			
		for (int i = 0; i < 3; i++) {
			stress_times_vel[i] = 0;
			for (int j = 0; j < 3; j++) {
				stress_times_vel[i] += stress_tens[i][j] * velF[j];
			}
			visc_stress_en_tot[i] = visc * stress_times_vel[i];
		}

		fluxViscEner[0] = visc_stress_en_tot[0]*nVec[0];
		fluxViscEner[1] = visc_stress_en_tot[1]*nVec[1];
		fluxViscEner[2] = visc_stress_en_tot[2]*nVec[2];
		
		areaF = area[ifc];	
		atomicAdd( &RHS[fc2el[ifc*2+0]*10+8], (fluxViscEner[0] + fluxViscEner[1] + fluxViscEner[2])*areaF);

		ifc += blockDim.x * gridDim.x;
	}
}
