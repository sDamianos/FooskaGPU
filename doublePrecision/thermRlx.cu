#include "strdata.h"
#include "eqos.h"

__device__ __noinline__ void initializeVariables(int *primMat, double avfMin, double *avf, double *rho, double *ar, double *enr, double *rhoMix, double *vel, double *Y, double *x, double *rE, double *abortX, double *consVec){
		
	avf[0] = fmin(fmax(avfMin/100,consVec[0]),1-avfMin/100);
	avf[1] = 1-avf[0];
	rho[0] = consVec[1]/avf[0];
	rho[1] = consVec[2]/avf[1];
	ar[0]  = avf[0]*rho[0];
	ar[1]  = avf[1]*rho[1];
	enr[0] = consVec[6]/ar[0];
	enr[1] = consVec[7]/ar[1];	
	*rhoMix = ar[0] + ar[1];
	vel[0] = consVec[3]/(*rhoMix);
	vel[1] = consVec[4]/(*rhoMix);
	vel[2] = consVec[5]/(*rhoMix);
	Y[0]   = ar[0]/(*rhoMix);
	Y[1]   = ar[1]/(*rhoMix);
	*rE    = consVec[8];
	
	x[0]   = rho[0];
	x[1]   = rho[1];
	x[2]   = enr[0];
	x[3]   = enr[1];

	abortX[0] = x[0];
	abortX[1] = x[1];
	abortX[2] = x[2];
	abortX[3] = x[3];
	
	double max=-1.0;
	for (int imat = 0; imat < 2; imat++) {
		if (avf[imat] > max) {
			max = avf[imat];
			*primMat = imat;
		}		
	}

} 

__device__ __noinline__ void fntPressRlx(INPUT *d_input, double avf[2], double rho[2], double ener[2], float *dP, float *dT){

	float de_drho[2], de_dP[2], p_I[2], c[2];

	float m = d_input->thermRlxRate;
	float press[2];
	press[0] = eqos(d_input,0,0,rho[0],ener[0]); 
	press[1] = eqos(d_input,1,0,rho[1],ener[1]); 
					        
	double kappa = 0;
	for(int iv = 0; iv < 2; iv++){
		de_drho[iv] = -( (d_input->gamma[iv]*d_input->pInf[iv] + press[iv])/((d_input->gamma[iv] - 1)*(rho[iv]*rho[iv])) );
		de_dP[iv]   = 1/((d_input->gamma[iv] - 1)*rho[iv]);
		p_I[iv] = press[iv];
		c[iv]   = ( p_I[iv]/(rho[iv]*rho[iv]) - de_drho[iv] )/ de_dP[iv];
		kappa  += ( rho[iv]*c[iv]/avf[iv] ); 
	}																											        
	kappa = kappa*(avf[1])*m;
	*dP = (press[0] - press[1])*exp(-fabs(kappa)*d_input->dt);
			
	float theta, temp[2];
	temp[0] = eqos(d_input,0,1,rho[0],ener[0]); 
	temp[1] = eqos(d_input,1,1,rho[1],ener[1]); 
	theta = (press[0]-press[1] > 1e-5) ? (*dP)/(press[0] - press[1]) : 0.0;
	*dT =  theta * (temp[0]-temp[1]);
}       


__device__ void computeResiduals(INPUT *d_input, int primMat, double x[4], double ar[2], double enr[2], float dP, float dT, double *res0, double *res1){

	double sum1,sum2;
	double press[2],temp[2];
	for(int imat = 0; imat < 2; imat++){
		if(imat!=primMat){
			sum1 = ar[imat]/x[imat];
			sum2 = ar[imat]*(x[2+imat] - enr[imat]);
		}
	}
	x[primMat]   = (1/(1-sum1)) * ar[primMat];
	x[2+primMat] = enr[primMat] - (1/ar[primMat])*sum2; 

	for(int imat = 0; imat < 2; imat++){
		press[imat] = eqos(d_input,imat,0,x[imat],x[imat+2]);
		temp[imat]  = eqos(d_input,imat,1,x[imat],x[imat+2]);
	}
	
	*res0 = press[0] - press[1] - dP;
	*res1 = temp[0]  - temp[1]  - dT;
}

__device__ __noinline__ void computeJacobian(INPUT *d_input, int primMat, double x[4], double ar[2], double enr[2], float dP, float dT, double res[2], double jac[2][2], int *abortFlag){
	
	double e[2];
	double res_df[2][2];
	int iv = 0;
	for(int j = 0; j < 4; j++){
		if(j!=primMat && j!=primMat+2){
			e[iv] = 1e-4*x[j];
			x[j] += e[iv];

			computeResiduals(d_input,primMat,x,ar,enr,dP,dT,&res_df[0][iv], &res_df[1][iv]);

			x[j] -= e[iv];
			iv++;
		}
	}

	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			jac[i][j] = (res_df[i][j]-res[i])/e[j];
			if(isnan(jac[i][j])==1){*abortFlag=1;}
		}
	}
}

__device__ __noinline__ void inverseJacobian(double jac[2][2], double invJac[2][2], int *abortFlag){

	double  det    = jac[0][0]*jac[1][1] - jac[1][0]*jac[0][1]; //maybe has to be a double
	invJac[0][0] = jac[1][1]/det;
	invJac[1][1] = jac[0][0]/det;
	invJac[1][0] = -jac[1][0]/det;
	invJac[0][1] = -jac[0][1]/det;
	if(det==0){*abortFlag=1;}

}

__device__ __noinline__ void abortCheckTherm(double *x, int abortFlag, double *abortX, float *minRes, double *res){

	if(abortFlag == 0 && x[0] > 0.0 && x[1] > 0.0){
		float check = sqrt(res[0]*res[0] + res[1]*res[1]);
		if(check < *minRes){
			*minRes = check;
			for(int i = 0; i < 4; i++) {
				abortX[i] = x[i];
			}
		}
	}

}

__device__ __noinline__ void updateSolution(bool *convergence, int itrMax, double urf, int primMat, int itr, double x[4], double res[2], double invJac[2][2], int abortFlag, double *abortX){

	int iv=0;
	double dx;
	for(int i = 0; i < 4; i++){
		if( i!=primMat && i!=primMat+2){
			dx = 0;
			for(int j = 0; j < 2; j++){
				dx += invJac[iv][j] * res[j];
			}
			x[i] = x[i] - urf*dx;
			iv++;
		}
	}

	int check = 0;
	for(int j = 0; j < 2; j++){
		if(fabs(res[j])>1.e-2){
			check += 1;
		}
	}
	if(check == 0){*convergence = true;} 
	if(abortFlag == 1 || itr > itrMax){
		for(int i = 0; i < 4; i++){
			x[i] = abortX[i];
		}
		*convergence = true;
	}
}

__device__ __noinline__ void updateConsVec(int primMat, double avfMin, double x[4], double ar[2], double enr[2], double *consVec){

	double sum1,sum2;
	for(int imat = 0; imat < 2; imat++){
		if(imat!=primMat){
			sum1 = ar[imat]/x[imat];
			sum2 = ar[imat]*(x[2+imat] - enr[imat]);
		}
	}
	x[primMat]   = (1/(1-sum1)) * ar[primMat];
	x[2+primMat] = enr[primMat] - (1/ar[primMat])*sum2; 

	consVec[0] = fmin(fmax(avfMin/100,ar[0]/x[0]),1-avfMin/100);
	consVec[6] = ar[0]*x[2];
	consVec[7] = ar[1]*x[3];
}

__global__ void thermalRlx(INPUT *d_input, int nElem, int itrMax, double urf,  int *idTherm2iel, double *consVec){
	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	int primMat,itr,iel;
	int abortFlag = 0;
	double avf[2],rho[2],ar[2],enr[2],vel[3],Y[2];
	double x[4],res[2],jac[2][2],invJac[2][2];
	double rhoMix,rE,de;
	double abortX[4];
	float minRes=1e9;
	bool convergence;

	while(tid < nElem){
		convergence = false;
		iel = idTherm2iel[tid];

		initializeVariables(&primMat,d_input->avfMin,avf,rho,ar,enr,&rhoMix,vel,Y,x,&rE,abortX,&consVec[iel*10]); 

		de = rE/rhoMix - Y[0]*enr[0] - Y[1]*enr[1] - 0.5*(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]); 
		enr[0] += de;
		enr[1] += de;

//		For inclusion of the finite relaxation methodology I have to compute pressure and include the dp in every residual computation, as it changes!!
		float dP = 0.0; float dT = 0.0;
		if(d_input->thermFntRlxFlag == 1){fntPressRlx(d_input,avf,rho,enr,&dP,&dT);}
		
		itr = 0;
		while(!convergence){
			itr++;	

			computeResiduals(d_input,primMat,x,ar,enr,dP,dT,&res[0],&res[1]); //include input struct!!!
		
			computeJacobian(d_input,primMat,x,ar,enr,dP,dT,res,jac,&abortFlag);

			inverseJacobian(jac,invJac,&abortFlag);

			abortCheckTherm(x,abortFlag,abortX,&minRes,res);

			updateSolution(&convergence,itrMax,urf,primMat,itr,x,res,invJac,abortFlag,abortX);
		}

		updateConsVec(primMat,d_input->avfMin,x,ar,enr,&consVec[iel*10]);

		tid += blockDim.x * gridDim.x;
	}
}
