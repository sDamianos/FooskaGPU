#include "strdata.h"
#include "eqos.h"

__device__ float saturationTemp(INPUT *d_input, float P, float temp, float urf){
	
	double fun_dt, grad_f, T;

	float Pinf 	= d_input->pInf[0];
	float g1   	= d_input->gamma[0]; 
	float g2   	= d_input->gamma[1];
	float Cp1  	= d_input->Cp[0];
	float Cp2  	= d_input->Cp[1];
	float Cv1  	= Cp1/g1;
	float Cv2  	= Cp2/g2;
	float eta1 	= d_input->eta[0]; 
	float eta2     = d_input->eta[1];
	float etaTone1 = 0;
	float etaTone2 = d_input->satConst[0]*log10(P)+d_input->satConst[1];
	float As       = (Cp1 - Cp2 + etaTone2 - etaTone1)/(Cp2 - Cv2);
	float Bs       = (eta1 - eta2)/(Cp2 - Cv2);
	float Cs       = (Cp2-Cp1)/(Cp2 - Cv2);
	float Ds       = (Cp1-Cv1)/(Cp2-Cv2);

	double dt = 1.e-3;
	double fun = 1;
	int itr = 0;
	
	T = temp;	
	while(fabs(fun)>1.e-2 && itr<50){
		fun = As +Bs/T + Cs*log10(T) + Ds*log10(P + Pinf) - log10(P);
		fun_dt = As + Bs/(T+dt) + Cs*log10(T+dt) + Ds*log10(P + Pinf) - log10(P);
		grad_f = (fun_dt - fun)/dt;
		T = T - urf*fun/grad_f;
		itr++;
	}

	return T;
}

// The cost of this functions is abou 4ms for 3 Million Cells. This consits the 
// 5 % of the total time step in first order RK or the 1.5% in 3rd order RK.
// The increament in the computation time of the relaxatio stpe is about 40%.
// However the whole cost of the function is to compute Tsat, and not the 
// atomicAdd. The Tsat was to be computed anyway and thus is a good practice
// to use this function as it almost cost 0 for building the arrays and are of good use.
__global__ void partitionRelaxation(INPUT *d_input, int nElem, float urf, int *idTherm2iel, int *idChem2iel, int *nTherm, int *nChem, double *consVec){

	int iel = threadIdx.x + blockIdx.x * blockDim.x;

	int idx;
	float  rho, enr, tempL, pressL, Tsat;	
	double avfMin = d_input->avfMin;
	double avf;

	while(iel < nElem){
		
		avf = fmin(fmax(avfMin/100,consVec[iel*10+0]),1-avfMin/100);
		rho = consVec[iel*10+1]/avf;
		enr = consVec[iel*10+6]/(avf*rho);

		tempL  = eqos(d_input,0,1,rho,enr);
		pressL = eqos(d_input,0,0,rho,enr);
		Tsat  = saturationTemp(d_input, pressL, tempL, urf);
//		if(iel == 29193){printf("Tsat = %lf\n", Tsat);}
	
		if(tempL > Tsat && avf > d_input->chemLimit){
			idx = atomicAdd(nChem, 1);
			idChem2iel[idx] = iel;
		}else{
			idx = atomicAdd(nTherm, 1);
			idTherm2iel[idx] = iel;
		}

		iel += gridDim.x * blockDim.x;
	}
}

void relaxation(MESH *mesh, FIELD *field, COMM *comm, cudaStream_t sMain, cudaStream_t sProx, cudaEvent_t evStepDone){
	
//	cudaEvent_t start, stop;
//	cudaEventCreate(&start);
//	cudaEventCreate(&stop);	
//	cudaEventRecord(start,sMain);	
	
	if(input.chemFlag == 1){
		cudaMemset(field->nTherm, 0, sizeof(int));
		cudaMemset(field->nChem, 0, sizeof(int));	
		partitionRelaxation<<<(mesh->nElem+511)/512,512,0,sMain>>>(d_input, mesh->nElem, input.satTempUrf, field->idTherm2iel, field->idChem2iel, field->nTherm, field->nChem, field->consVec); 
	}
	
	cudaStreamSynchronize(sMain);
	
	int h_nTherm, h_nChem;
	cudaMemcpy(&h_nTherm, field->nTherm, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&h_nChem,  field->nChem,  sizeof(int), cudaMemcpyDeviceToHost);
	
//	printf("nElem nThermCells nChemCells %d %d %d\n", mesh->nElem, h_nTherm, h_nChem);
	
	thermalRlx<<<(h_nTherm+511)/512,512,0,sMain>>>(d_input,h_nTherm,input.thermItrMax,input.thermUrf,field->idTherm2iel,field->consVec);

	if(input.chemFlag==1){
		int nProc = max(h_nChem,1);
		chemRlx<<<(nProc+511)/512,512,0,sMain>>>(d_input,h_nChem,input.chemItrMax,input.chemUrf,field->idChem2iel,field->consVec);
	}
	
//	cudaEventRecord(stop, sMain);
//	cudaEventSynchronize(stop);
//	float   elapsedTime;
//	cudaEventElapsedTime(&elapsedTime,start,stop);
//	printf("relaxation time %3.5f ms\n", elapsedTime);
	
	cudaEventRecord(evStepDone, sMain);
	cudaStreamWaitEvent(sProx, evStepDone, 0);
	
}
