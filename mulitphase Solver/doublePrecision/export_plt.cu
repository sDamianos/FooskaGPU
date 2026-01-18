#include "strdata.h"
#include "eqos.h"

__constant__ int TPNDS[8][9] = {
    {0, 0, 0, 0, 1, 0, 2, 0, 0},   
    {0, 0, 0, 0, 2, 1, 5, 0, 1},  
    {0, 0, 0, 0, 2, 3, 3, 0, 3},   
    {0, 0, 0, 0, 0, 2, 0, 0, 2},  
    {0, 0, 0, 0, 3, 4, 1, 0, 4},   
    {0, 0, 0, 0, 3, 4, 4, 0, 5},   
    {0, 0, 0, 0, 3, 4, 4, 0, 7},  
    {0, 0, 0, 0, 3, 4, 1, 0, 6}    
};

__global__ void computeConnectivityVector(int rank, int nElem, int *lcl2glbEl, int *elNd2lclNd, int *lcl2glbNd, int *nElNds, int *connecVec){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int iel = tid/8;
	int ind = tid%8;

	int lclNd, glbNd, ielGlb;
	
	while(iel < nElem){
		lclNd  = elNd2lclNd[iel*8 + TPNDS[ind][nElNds[iel]]];
		glbNd  = lcl2glbNd[lclNd];
		ielGlb = lcl2glbEl[iel]; 
		connecVec[ielGlb*8 + ind] = 1 + glbNd; 
//		if(rank==2){printf("iel %d ind %d lclNd %d glbNd %d ielGlb %d\n", iel, ind, lclNd, glbNd, ielGlb);}
	
		tid += blockDim.x * gridDim.x;
		iel  = tid/8;
		ind  = tid%8;
	}
}
		

__global__ void computeNodeVector(INPUT *d_input,int rank, int nNodes, int maxEl4nd, int nOut, int *lcl2glbNd, int *nEl4nd, int *idEl4nd, double *consVec, double *grad, double *glbNdCord, double *vecNd){

	int ilnd = threadIdx.x + blockIdx.x * blockDim.x;

	double avfL,avfG,rhoL,rhoG,rho;
	double u,v,w,eIntL,eIntG;
	double press,temp,vapM;
	double gradRho,dRhoDx,dRhoDy,dRhoDz;
	int iel,ignd;

	while(ilnd < nNodes){
		for(int iv = 0; iv < nEl4nd[ilnd]; iv++){
			iel = idEl4nd[ilnd*maxEl4nd + iv];
			if(iel > -1 ){
				
				avfL  = consVec[iel*d_input->NEQ+0];
				avfG  = 1 - avfL;
				rhoL  = consVec[iel*d_input->NEQ+1]/avfL;
				rhoG  = consVec[iel*d_input->NEQ+2]/avfG;
				rho   = avfL*rhoL + avfG*rhoG;
				
				u     = consVec[iel*d_input->NEQ+3]/rho;
				v     = consVec[iel*d_input->NEQ+4]/rho;
				w     = consVec[iel*d_input->NEQ+5]/rho;
				
				eIntL  = consVec[iel*d_input->NEQ+6]/(avfL*rhoL);
				eIntG  = consVec[iel*d_input->NEQ+7]/(avfG*rhoG);
			
				press  = avfL * eqos(d_input,0,0,rhoL,eIntL) + avfG * eqos(d_input,1,0,rhoG,eIntG);
				temp   = avfL * eqos(d_input,0,1,rhoL,eIntL) + avfG * eqos(d_input,1,1,rhoG,eIntG);
				vapM   = consVec[iel*d_input->NEQ+9];

				dRhoDx = grad[iel*d_input->NEQ*3 + 9*3 + 0];
				dRhoDy = grad[iel*d_input->NEQ*3 + 9*3 + 1];
				dRhoDz = grad[iel*d_input->NEQ*3 + 9*3 + 2];
				gradRho = sqrt(dRhoDx*dRhoDx + dRhoDy*dRhoDy + dRhoDz*dRhoDz);
				
				ignd = lcl2glbNd[ilnd];
				vecNd[ignd*nOut + 0 ] += glbNdCord[ilnd*3 + 0]/nEl4nd[ilnd]; //x
				vecNd[ignd*nOut + 1 ] += glbNdCord[ilnd*3 + 1]/nEl4nd[ilnd]; //y
				vecNd[ignd*nOut + 2 ] += glbNdCord[ilnd*3 + 2]/nEl4nd[ilnd]; //z
	
				vecNd[ignd*nOut + 3 ] += avfL/nEl4nd[ilnd];
				vecNd[ignd*nOut + 4 ] += avfG/nEl4nd[ilnd];
				vecNd[ignd*nOut + 5 ] += rho/nEl4nd[ilnd];
				vecNd[ignd*nOut + 6 ] += u/nEl4nd[ilnd];
				vecNd[ignd*nOut + 7 ] += v/nEl4nd[ilnd];
				vecNd[ignd*nOut + 8 ] += w/nEl4nd[ilnd];
				vecNd[ignd*nOut + 9 ] += press/nEl4nd[ilnd];	
				vecNd[ignd*nOut + 10] += temp/nEl4nd[ilnd];	
				vecNd[ignd*nOut + 11] += vapM/nEl4nd[ilnd];	
				vecNd[ignd*nOut + 12] += gradRho/nEl4nd[ilnd];	
			}
		}
		ilnd += blockDim.x * gridDim.x;
	}
}

void export_plt(EXPORT *exp, FIELD *field){
	
	int nOut,rank;
	int * h_connecVec, *connecVec;
	double * h_vecNd, *vecNd;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	int nNodes = exp->nNodes;
	int nElem = exp->nElem;
	int nGlbNodes = exp->nGlbNodes;
	int nGlbElem  = exp->nGlbElem;
	
	nOut = 13;
	
	h_vecNd = (double*) malloc(nGlbNodes*nOut*sizeof(double));
	h_connecVec = (int*) malloc(nGlbElem*8*sizeof(int));
	
	for(int i = 0; i < nGlbNodes*nOut; i++){
		h_vecNd[i] = 0.0;
	}
	for(int i = 0; i < nGlbElem*8; i++){
		h_connecVec[i] = 0;
	}
	
	cudaMalloc((void**)&vecNd, nGlbNodes  * nOut * sizeof(double));
	cudaMalloc((void**)&connecVec, nGlbElem * 8 * sizeof(int));
		
	cudaMemcpy(vecNd,     h_vecNd,      nGlbNodes * nOut * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(connecVec, h_connecVec,  nGlbElem * 8     * sizeof(int),   cudaMemcpyHostToDevice);
		
	computeNodeVector<<<(nNodes*nOut+511)/512,512>>>(d_input,rank,nNodes,exp->maxEl4nd,nOut,exp->lcl2glbNd,exp->nEl4nd,exp->idEl4nd,field->consVec,field->grad,exp->glbNdCord,vecNd);
	computeConnectivityVector<<<(nElem*8+511)/512,512>>>(rank,nElem,exp->lcl2glbEl,exp->elNd2lclNd,exp->lcl2glbNd,exp->nElNds,connecVec);
	
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("CUDA error in export: %s\n", cudaGetErrorString(err));
	}

	cudaDeviceSynchronize();
	MPI_Allreduce(MPI_IN_PLACE, vecNd, nGlbNodes*nOut, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, connecVec, nGlbElem*8, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			
	cudaMemcpy(h_vecNd,      vecNd,     nGlbNodes * nOut * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_connecVec,  connecVec, nGlbElem  * 8    * sizeof(int),   cudaMemcpyDeviceToHost);
	
	cudaDeviceSynchronize();
	MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0){
  		int i,j;
  		int len = 0;
		FILE *f1;
        	char file[50];
        	sprintf(file,"field_%07d.dat",iStep+1);
        	f1 = fopen(file,"w");
		fprintf(f1, "VARIABLES = \"X\" \"Y\" \"Z\" \"avfL\" \"avfG\" \"rho\" \"U\" \"V\" \"W\" \"Pressure\" \"Temperature\" \"vaporMass\" \"gradRho\"\n");
		fprintf(f1, "ZONE T=\"1\",N=%d,E=%d,F=FEPOINT,ET=BRICK,DT=(double,double,double),SOLUTIONTIME=%le\n",  nGlbNodes, nGlbElem, physicalTime);
		for(int ind = 0; ind < nGlbNodes; ind++){
			for(int i = 0; i < nOut; i++){
				fprintf(f1, "%le ", h_vecNd[ind*nOut + i]);
			}
			 fprintf(f1, "\n");
		}
		for(int iel = 0; iel < nGlbElem; iel++){
			for(int ilnd = 0; ilnd < 8; ilnd++){
				fprintf(f1, "%d ", h_connecVec[iel*8 + ilnd]);
			}
			 fprintf(f1, "\n");
		}
		fclose(f1);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	free(h_vecNd);
	free(h_connecVec);
	cudaFree(connecVec);
	cudaFree(vecNd);
	
}
/*
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("CUDA error: %s\n", cudaGetErrorString(err));
	}
*/
