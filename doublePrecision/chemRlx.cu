#include "strdata.h"
#include "eqos.h"

__device__ __noinline__ void initializeVariables(int *primMat, double avfMin, double *rhoMix, double *rEint, double *x, double *abortX, double *consVec){
		
	double u,v,w;

	x[0] = fmin(fmax(avfMin/100,consVec[0]),1-avfMin/100);
	x[1] = 1-x[0];
	x[2] = consVec[1]/x[0]; 
	x[3] = consVec[2]/x[1];
	x[4] = consVec[6]/(x[0]*x[2]);
	x[5] = consVec[7]/(x[1]*x[3]);

	abortX[0] = x[0];
	abortX[1] = x[1];
	abortX[2] = x[2];
	abortX[3] = x[3];
	abortX[4] = x[4];
	abortX[5] = x[5];
	
	*rhoMix = x[0]*x[2] + x[1]*x[3];
	u = consVec[3]/(*rhoMix);
	v = consVec[4]/(*rhoMix);
	w = consVec[5]/(*rhoMix);
	*rEint = consVec[8] - 0.5*(*rhoMix)*(u*u+v*v+w*w);

	double max=-1.0;
	for (int imat = 0; imat < 2; imat++) {
		if (x[imat] > max) {
			max = x[imat];
			*primMat = imat;
		}		
	}

} 

__device__ double gibbsEnr(INPUT *d_input, int imat, double val1, double val2){

	double eta_tone[2];

	eta_tone[0] = 0;
	eta_tone[1] = d_input->satConst[0]*log10(val1)+d_input->satConst[1];
	
	double enthalpy = d_input->Cp[imat]*val2+d_input->eta[imat];
	double entropy = (d_input->Cp[imat]/d_input->gamma[imat])*log10(pow(val2, d_input->gamma[imat])/pow(val1 + d_input->pInf[imat], d_input->gamma[imat]-1))+eta_tone[imat];

	double gibbs = enthalpy - val2*entropy;
	return gibbs;
}

__device__ void computeResiduals(INPUT *d_input, int primMat, double x[6], double rhoMix, double rEint, double *res0, double *res1, double *res2){

	double sum1, sum2, sum3;
	double press[2], temp[2], gibbs[2];
	for(int imat = 0; imat < 2; imat++){
		if(imat!=primMat){
			sum1 = x[imat];
			sum2 = x[imat]*x[imat+2];    
			sum3 = x[imat]*x[imat+2]*x[imat+4];
		}
	}
	x[primMat]   = 1 - sum1;
	x[primMat+2] = (rhoMix - sum2)/x[primMat];
	x[primMat+4] = (rEint - sum3)/(x[primMat]*x[primMat+2]);

	for(int imat = 0; imat < 2; imat++){
		press[imat] = eqos(d_input,imat,0,x[imat+2],x[imat+4]); 
		temp[imat]  = eqos(d_input,imat,1,x[imat+2],x[imat+4]);
		gibbs[imat] = gibbsEnr(d_input, imat, press[imat], temp[imat]); 
	}

	*res0 = press[0] - press[1];
	*res1 = temp[0]  - temp[1];
	*res2 = gibbs[0] - gibbs[1];	
}

__device__ __noinline__ void computeJacobian(INPUT *d_input, int primMat, double x[6], double rhoMix, double rEint, double res[3], double jac[3][3], int *abortFlag){

	double e[3];
	double res_df[3][3];
	int iv = 0;
	for(int j = 0; j < 6; j++){
		if(j!=primMat && j!=primMat+2 && j!=primMat+4){
			e[iv] = fabs(x[j]) < 1e-9 ? 1e-10 : 1.e-5*x[j];
			x[j] += e[iv];
			
			computeResiduals(d_input,primMat,x,rhoMix,rEint,&res_df[0][iv],&res_df[1][iv],&res_df[2][iv]); 
			
			x[j] -= e[iv];
			iv++;
		}
	}

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			jac[i][j] = (res_df[i][j]-res[i])/e[j];
			if(isnan(jac[i][j])==1){*abortFlag=1;}
		}
	}
}

__device__ __noinline__ void inverseJacobian(double jac[3][3], double invJac[3][3],int *abortFlag){

	double det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
		   - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
		   + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);

	if(det == 0){*abortFlag = 1;}

	double inv_det = 1.0 / det;
	invJac[0][0] =  inv_det * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]);
	invJac[0][1] = -inv_det * (jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1]);
	invJac[0][2] =  inv_det * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);

	invJac[1][0] = -inv_det * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]);
	invJac[1][1] =  inv_det * (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]);
	invJac[1][2] = -inv_det * (jac[0][0] * jac[1][2] - jac[0][2] * jac[1][0]);

	invJac[2][0] =  inv_det * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
	invJac[2][1] = -inv_det * (jac[0][0] * jac[2][1] - jac[0][1] * jac[2][0]);
	invJac[2][2] =  inv_det * (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]);
}

__device__ __noinline__ void abortCheck(double *x, int abortFlag, double *abortX, float *minRes, double *res){

	if(abortFlag == 0 && x[0]>0.0 && x[0]<1.0 && x[1]>0.0 && x[1] < 1.0 && x[2] > 0.0 && x[3] > 0.0){ //for laser ablation extra checks are needed!!!
		float check = sqrt(res[0]*res[0]+res[1]*res[1]+res[2]*res[2]);
		if( check < *minRes){
			*minRes = check;
			for(int i = 0; i < 6; i++) {
				abortX[i] = x[i];
			}
		}
	}	

}

__device__ __noinline__ void checkUpdate(INPUT *d_input, double xOld, double *x, double dx, int i, float urf, int *abortFlag, double rho, double rEint){

	if(i==0||i==1){
		if(x[i]<0){while(x[i] < d_input->chemLimit     && urf > 0.0005){urf = urf/2; x[i] = xOld - urf*dx;}}
		if(x[i]>1){while(x[i] > 1.0-d_input->chemLimit && urf > 0.0005){urf = urf/2; x[i] = xOld - urf*dx;}}
		if(x[i]<0){*abortFlag = 1;}
		if(x[i]>1){*abortFlag = 1;}
	}
	else if(i==2||i==3){
		double ar = x[i-2]*x[i];
		if(x[i] < 0){  while(x[i] < 0   && urf > 0.0005){urf = urf/2; x[i] = xOld - urf*dx;}}
		if(ar   > rho){while(ar   > rho && urf > 0.0005){urf = urf/2; x[i] = xOld - urf*dx; ar = x[i-2]*x[i];}}
		if(x[i] < 0){  *abortFlag = 1;}
		if(ar   > rho){*abortFlag = 1;}
	}
	else{
	}
//	else{  //if enery is negative can  not stand
//		double are = x[i-4]*x[i-2]*x[i];
//		if(are > rEint){while(are > rEint && urf > 0.0005){urf = urf/2; x[i] = xOld - urf*dx; are = x[i-4]*x[i-2]*x[i];}}
//		if(are > rEint){*abortFlag =1;}
//	}

}

__device__ __noinline__ void updateSolution(INPUT *d_input, bool *convergence, int itrMax, double urf, int primMat, int itr, double x[6], double res[3], double invJac[3][3], int abortFlag, double *abortX, double rhoMix, double rEint){

	int iv = 0;
	double dx;
	for(int i = 0; i < 6; i++){
		if(i!=primMat && i!=primMat+2 && i!=primMat+4){
			dx = 0;
			for (int j = 0; j < 3; j++) {
				dx += invJac[iv][j] * res[j];
			}
			x[i] = x[i] - urf*dx;
			checkUpdate(d_input,x[i]+urf*dx,x,dx,i,urf,&abortFlag,rhoMix,rEint);
			iv++;
		}
	}
	
	int check = 0;
	for(int j = 0; j < 3; j++){
		if(fabs(res[j])>1.e-1){
			check += 1;
		}
	}
	if(check == 0){*convergence = true;} 
	if(abortFlag == 1 || itr > itrMax){
		for(int i = 0; i < 6; i++){
			x[i] = abortX[i];
		}
		*convergence = true;
	}
}

__device__ __noinline__ void updateConsVec(int primMat, double avfMin, double x[6], double rhoMix, double rEint, double *consVec){

	double sum1, sum2, sum3;
	for(int imat = 0; imat < 2; imat++){
		if(imat!=primMat){
			sum1 = x[imat];
			sum2 = x[imat]*x[imat+2];    
			sum3 = x[imat]*x[imat+2]*x[imat+4];
		}
	}
	x[primMat]   = 1 - sum1;
	x[primMat+2] = (rhoMix - sum2)/x[primMat];
	x[primMat+4] = (rEint - sum3)/(x[primMat]*x[primMat+2]);

	consVec[0] = fmin(fmax(avfMin/100,x[0]),1-avfMin/100);
	consVec[1] = x[0]*x[2];
	consVec[2] = x[1]*x[3];
	consVec[6] = x[0]*x[2]*x[4];
	consVec[7] = x[1]*x[3]*x[5];
	consVec[9] = x[1]*x[3];
}

__global__ void chemRlx(INPUT *d_input, int nElem, int itrMax, double urf, int *idChem2iel, double *consVec){
	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	int primMat,itr,iel;
	int abortFlag = 0;
	double x[6],res[3],jac[3][3],invJac[3][3];
	double rhoMix,rEint;
	double abortX[6];
	float minRes=1e9;
	bool convergence;

	while(tid < nElem){
		convergence = false;
		iel = idChem2iel[tid];

		initializeVariables(&primMat,d_input->avfMin,&rhoMix,&rEint,x,abortX,&consVec[iel*10]);
		// re_inti is needed only in case of finite mass transfer to compute g1 - g2
		// with the maximum accuracy and consistency to total energy. Ohterwise energy
		// equilibrium is guranteed by the relaxation step itself.

//		finite mass transfer functions!!
//		dG = (1-d_input->fntMassTransfVel
//		dP = fntPressRelax;
//		dT = as in thermalRlx		

		//if(tid==0){printf("avf %le %le rho %f %f enr %f %f arg %f abort %d\n", x[0], x[1], x[2], x[3], x[4], x[5], consVec[iel*10+9], abortFlag);}
		
		itr = 0;
		while(!convergence){
			itr++;	

			computeResiduals(d_input,primMat,x,rhoMix,rEint,&res[0],&res[1],&res[2]); 
			
			computeJacobian(d_input,primMat,x,rhoMix,rEint,res,jac,&abortFlag);

			inverseJacobian(jac,invJac,&abortFlag);

			abortCheck(x,abortFlag,abortX,&minRes,res);
			
//			if(tid==0){printf("itr %d avf %e %e rho %f %f enr %f %f res %e %e %e abort %d\n", itr,  x[0], x[1], x[2], x[3], x[4], x[5], res[0], res[1], res[2], abortFlag);}

			updateSolution(d_input,&convergence,itrMax,urf,primMat,itr,x,res,invJac,abortFlag,abortX,rhoMix,rEint);
		}

		updateConsVec(primMat,d_input->avfMin,x,rhoMix,rEint,&consVec[iel*10]);
	//	if(tid==0){printf("Final: avf %le %le rho %f %f enr %f %f arg %f abort %d\n", x[0], x[1], x[2], x[3], x[4], x[5], consVec[iel*10+9], abortFlag);}

		tid += blockDim.x * gridDim.x;
	}
}
