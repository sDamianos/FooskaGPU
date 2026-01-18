#include "strdata.h"

void partitionConnectivity(int nLclElem, int rank, int size, int nNeigRanks, int max, int *nFaces, int *lcl2glbEl, int *lcl2glbFc, int *glb2lclFc, int *neig, int *part, int *elFc2fc, int *neigRanks, int *nProxFaces, int *neigRank4fc, int *lclProx4fc, int *lclFc2idRcv){
	
	int *sendBuff_idProx2glbFc = (int*)malloc(size*max*sizeof(int)); 
	int *recvBuff_idProx2glbFc = (int*)malloc(size*max*sizeof(int));
	int *iProx                 = (int*)malloc(size*sizeof(int));
	
	for(int i = 0; i < size; i++){
		iProx[i] = 0;
	}

	int ielGlb, ifcGlb, irank, neigRank;
	int glbFc, lclFc, ifcLcl;
	
	for(int iel = 0; iel < nLclElem; iel++){
		ielGlb = lcl2glbEl[iel];
		for(int ifc = 0; ifc < nFaces[ielGlb]; ifc++){
			ifcLcl = elFc2fc[iel*6+ifc];
			ifcGlb = lcl2glbFc[ifcLcl];
			if(neig[6*ielGlb + ifc] >= 0){
				neigRank = part[neig[6*ielGlb + ifc]];
				if(neigRank != rank){
					neigRank4fc[ifcLcl]   = neigRank;
					lclProx4fc[ifcLcl]   = iProx[neigRank];
					sendBuff_idProx2glbFc[neigRank*max + iProx[neigRank]] = ifcGlb;
					iProx[neigRank]++;
				}
			}
		}
	}

	MPI_Request reqs[2*nNeigRanks];
	int r=0;
	for(int ipr=0; ipr<nNeigRanks; ipr++){
		neigRank = neigRanks[ipr];
		MPI_Isend(&sendBuff_idProx2glbFc[neigRank*max], nProxFaces[neigRank], MPI_INT, neigRank, rank, MPI_COMM_WORLD, &reqs[r++]);
	}

	// Post all receives
	for(int ipr=0; ipr<nNeigRanks; ipr++){
		neigRank = neigRanks[ipr];
		MPI_Irecv(&recvBuff_idProx2glbFc[neigRank*max], nProxFaces[neigRank], MPI_INT, neigRank, neigRank, MPI_COMM_WORLD, &reqs[r++]);
	}

	MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);

	for(int ipr = 0; ipr < nNeigRanks; ipr++){
		irank = neigRanks[ipr];
		for(int iv = 0; iv < nProxFaces[irank]; iv++){
			glbFc = recvBuff_idProx2glbFc[irank*max + iv];
			lclFc = glb2lclFc[glbFc];
			lclFc2idRcv[lclFc] = iv;
		}
	}
	
	free(sendBuff_idProx2glbFc);
	free(recvBuff_idProx2glbFc);
	free(iProx);
}
