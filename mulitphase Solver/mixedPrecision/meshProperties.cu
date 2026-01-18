#include "strdata.h"

int h_nFaceNodes(int ifc, int nLclNodes){

  if      (nLclNodes == 8) 	      { return 4; }
  else if (nLclNodes == 4)	      { return 3; }
  else if (nLclNodes == 6 && ifc < 3) { return 4; }
  else if (nLclNodes == 6 && ifc > 2) { return 3; }
  else if (nLclNodes == 5 && ifc < 1) { return 4; }
  else if (nLclNodes == 5 && ifc > 0) { return 3; }
  else                 		      { return 0; }
}


int fcNd2ElNd(int ifc, int ifcnd, int nLclNodes){ 

int M[6][4];  
  if (nLclNodes==8) {
    
    M[0][0]=0; M[1][0]=3; M[2][0]=2; M[3][0]=2; M[4][0]=2; M[5][0]=6;
    M[0][1]=1; M[1][1]=7; M[2][1]=6; M[3][1]=0; M[4][1]=3; M[5][1]=4;
    M[0][2]=5; M[1][2]=5; M[2][2]=7; M[3][2]=4; M[4][2]=1; M[5][2]=5;
    M[0][3]=4; M[1][3]=1; M[2][3]=3; M[3][3]=6; M[4][3]=0; M[5][3]=7;
  
  } else if (nLclNodes==4) {

    M[0][0]=1; M[1][0]=0; M[2][0]=1; M[3][0]=2; M[4][0]=-1; M[5][0]=-1;
    M[0][1]=0; M[1][1]=1; M[2][1]=2; M[3][1]=0; M[4][1]=-1; M[5][1]=-1;
    M[0][2]=2; M[1][2]=3; M[2][2]=3; M[3][2]=3; M[4][2]=-1; M[5][2]=-1;
    M[0][3]=4; M[1][3]=5; M[2][3]=6; M[3][3]=7; M[4][3]=-1; M[5][3]=-1;

  } else if (nLclNodes==6){
    
    M[0][0]=0; M[1][0]=2; M[2][0]=0; M[3][0]=0; M[4][0]=3; M[5][0]=-1;
    M[0][1]=1; M[1][1]=5; M[2][1]=3; M[3][1]=2; M[4][1]=4; M[5][1]=-1;
    M[0][2]=4; M[1][2]=4; M[2][2]=5; M[3][2]=1; M[4][2]=5; M[5][2]=-1;
    M[0][3]=3; M[1][3]=1; M[2][3]=2; M[3][3]=6; M[4][3]=7; M[5][3]=-1;
  
  } else if (nLclNodes==5){
  
    
    M[0][0]=0; M[1][0] = 0; M[2][0] = 1; M[3][0] = 3; M[4][0] = 2; M[5][0] =-1;
    M[0][1]=2; M[1][1] = 1; M[2][1] = 3; M[3][1] = 2; M[4][1] = 0; M[5][1] =-1;
    M[0][2]=3; M[1][2] = 4; M[2][2] = 4; M[3][2] = 4; M[4][2] = 4; M[5][2] =-1;
    M[0][3]=1; M[1][3] =-1; M[2][3] =-1; M[3][3] =-1; M[4][3] =-1; M[5][3] =-1;
      
  }
  return M[ifc][ifcnd];
}

void normalVector(float *x, float *y, float *z, float *n){

	double norm;

	norm = sqrt(pow(( (y[1] - y[0]) * (z[2] - z[0]) - (y[2] - y[0]) * (z[1] - z[0])),2.0) + 
		    pow(( (x[1] - x[0]) * (z[2] - z[0]) - (x[2] - x[0]) * (z[1] - z[0])),2.0) + 
		    pow(( (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])),2.0));
    
	n[0] =  ((y[1] - y[0]) * (z[2] - z[0]) - (y[2] - y[0]) * (z[1] - z[0]))/norm;
	n[1] = -((x[1] - x[0]) * (z[2] - z[0]) - (x[2] - x[0]) * (z[1] - z[0]))/norm;
	n[2] =  ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))/norm;

}

void tangentVector(float *x, float *y, float *z, float *n, float *nt1, float *nt2){

	double norm, norm1;
	
	norm = sqrt( pow(x[1] - x[0],2.0) + pow(y[1] - y[0],2.0) + pow(z[1] - z[0],2.0) ); 

	nt1[0] = (x[1] - x[0])/norm;
	nt1[1] = (y[1] - y[0])/norm;
	nt1[2] = (z[1] - z[0])/norm;

	norm1 = sqrt( pow(n[1]*nt1[2] - nt1[1]*n[2],2.0) + 
		      pow(n[0]*nt1[2] - nt1[0]*n[2],2.0) + 
		      pow(n[0]*nt1[1] - nt1[0]*n[1],2.0) );

	nt2[0] =  (n[1]*nt1[2] - nt1[1]*n[2])/norm1;
	nt2[1] = -(n[0]*nt1[2] - nt1[0]*n[2])/norm1;
	nt2[2] =  (n[0]*nt1[1] - nt1[0]*n[1])/norm1;
}
			
void computeArea(int nFcNds, float *x, float *y, float *z, float *area){			

	double dx1,dy1,dz1;
	double dx2,dy2,dz2;
	
	area[0] = 0.0;
	for(int it = 0; it < nFcNds-2; it++){  // Loop for triangles
		dx1 = x[1+it] - x[0];
		dx2 = x[2+it] - x[0];
		dy1 = y[1+it] - y[0];
		dy2 = y[2+it] - y[0];
		dz1 = z[1+it] - z[0];
		dz2 = z[2+it] - z[0];
		area[0]+=fabs(0.5*sqrt(pow(dy1*dz2-dz1*dy2,2)+pow(dz1*dx2-dx1*dz2,2)+pow(dx1*dy2-dy1*dx2,2)));
	}
	
}
			
void computeVolume(int nElem, int iel, const int (*vnHexa)[3], const int (*vnTetra)[3], const int (*vnPrism)[3], const int (*vnPyramid)[3], int *nLclNodes, int *elNd2lclNd, float *glbNdCord, float *volume){
	

	int num_tetra;
	const int (*vn)[3];
	double dx[3],dy[3],dz[3];
	double cross[3];

	if (nLclNodes[iel] == 8){
		num_tetra=6; 
		vn = vnHexa;
	} else if (nLclNodes[iel] == 4){
		num_tetra=1;
		vn = vnTetra;
	} else if (nLclNodes[iel] == 6){
		num_tetra=3;
		vn = vnPrism;
	} else if (nLclNodes[iel] == 5){
		num_tetra=2;
		vn = vnPyramid;
	}

	volume[iel]=0.0;
	for(int it=0; it<num_tetra;it++){  
		for(int in=0; in<3;in++){  
			dx[in] = glbNdCord[elNd2lclNd[iel*8 + vn[it][in]]*3 + 0] - glbNdCord[elNd2lclNd[iel*8 + 0]*3 + 0];
			dy[in] = glbNdCord[elNd2lclNd[iel*8 + vn[it][in]]*3 + 1] - glbNdCord[elNd2lclNd[iel*8 + 0]*3 + 1];
			dz[in] = glbNdCord[elNd2lclNd[iel*8 + vn[it][in]]*3 + 2] - glbNdCord[elNd2lclNd[iel*8 + 0]*3 + 2];
		}	
		 
		cross[0]=dy[0]*dz[1]-dz[0]*dy[1];
		cross[1]=dz[0]*dx[1]-dx[0]*dz[1];
		cross[2]=dx[0]*dy[1]-dy[0]*dx[1];
		volume[iel]+=(cross[0]*dx[2]+cross[1]*dy[2]+cross[2]*dz[2])/6.0;    
	}
}

void meshProperties(HOSTLOCALDATA *hLclData, int *nLclNodes, int *fc2el, int *fc2elFc, int *elNd2lclNd){
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int iel, elFc, elNd, ind, nFcNds;
	float x[4],y[4],z[4];

	hLclData->n    = (float*)malloc(hLclData->nFcs*3*sizeof(float));
	hLclData->nt1  = (float*)malloc(hLclData->nFcs*3*sizeof(float));
	hLclData->nt2  = (float*)malloc(hLclData->nFcs*3*sizeof(float));
	hLclData->area = (float*)malloc(hLclData->nFcs*  sizeof(float));
	for(int ifc = 0; ifc < hLclData->nFcs; ifc++){
	
		iel = fc2el[ifc*2+0]; 
		elFc = fc2elFc[ifc*2+0];
		nFcNds = h_nFaceNodes(elFc,nLclNodes[iel]);
		for(int iFcNd = 0; iFcNd < nFcNds; iFcNd++){
			elNd = fcNd2ElNd(elFc,iFcNd,nLclNodes[iel]);
			ind  = elNd2lclNd[iel*8 + elNd];
			x[iFcNd] = hLclData->glbNdCord[ind*3 + 0];
			y[iFcNd] = hLclData->glbNdCord[ind*3 + 1];  
			z[iFcNd] = hLclData->glbNdCord[ind*3 + 2];  
		}
		
		normalVector (x, y, z, &hLclData->n[ifc*3]);
		tangentVector(x, y, z, &hLclData->n[ifc*3], &hLclData->nt1[ifc*3], &hLclData->nt2[ifc*3]);
		computeArea  (nFcNds, x, y, z, &hLclData->area[ifc]);
	}

	const int vnHexa[6][3] = { {1,3,5}, {5,3,7}, {5,7,4}, {3,2,7}, {7,2,6}, {7,6,4} };
	const int vnTetra[1][3] = { {1,2,3} };
	const int vnPrism[3][3] = { {1,2,5}, {1,5,4}, {4,5,3} };
	const int vnPyramid[2][3] = { {1,3,4}, {3,2,4} };
	
	hLclData->volume = (float*)malloc((hLclData->nElem+hLclData->nGhostElem)*sizeof(float));
	for(int iel = 0; iel < hLclData->nElem; iel++){
		computeVolume(hLclData->nElem, iel, vnHexa, vnTetra, vnPrism, vnPyramid, nLclNodes, elNd2lclNd, hLclData->glbNdCord, hLclData->volume);
	}
	
}
