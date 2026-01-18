#include "strdata.h"

int countLocalElem(int nGlbElem, int rank, int *part){
	
	int count = 0;
	for(int iel = 0; iel < nGlbElem; iel++){
		if(part[iel] == rank){
			count++;
		}	
	}
	return count;
}
	
int countGhostCells(int nGlbElem, int nLclElem, int rank, int *nOuter, int *nFaces, int *part, int *neig, int *outerFlag4glbEl){
	
	int iel1;
	int nGhostCells;
	int flag;
	int *seenGhostCell = (int*)malloc(nGlbElem*sizeof(int));
	for(int ielGlb = 0; ielGlb < nGlbElem; ielGlb++){
		seenGhostCell[ielGlb] = 0;
	}

	nGhostCells = 0;
	*nOuter = 0;
	for(int iel = 0; iel < nGlbElem; iel++){
		flag = 0;
		outerFlag4glbEl[iel] = 0;
		if(part[iel] == rank){
			for(int ifc = 0; ifc < nFaces[iel]; ifc++){
				iel1 = neig[iel*6 + ifc];
				if(iel1 != -1 && part[iel1] != rank && seenGhostCell[iel1] == 0){
					seenGhostCell[iel1] = 1;
					nGhostCells++;
				}
				if(iel1 != -1 && part[iel1] != rank){
					flag = 1;
				}
			}
			if(flag == 1){
				(*nOuter)++;
				outerFlag4glbEl[iel] = 1;
			}
		}
	}
	free(seenGhostCell);

	return nGhostCells;
}
	
void defElemPartitioning(int nElemGlb, int nLclElem, int nOuter, int rank, int *part, int *lcl2glbEl, int *glb2lclEl, int *outerFlag4glbEl){
	
	int ielLcl = 0;
	int ielOuter = nLclElem - nOuter;
	for(int ielGlb = 0; ielGlb < nElemGlb; ielGlb++){
		glb2lclEl[ielGlb] = -1;
		if(part[ielGlb] == rank && outerFlag4glbEl[ielGlb] == 0){
			lcl2glbEl[ielLcl] = ielGlb;
			glb2lclEl[ielGlb]  = ielLcl;
			ielLcl++;
		}else if(part[ielGlb] == rank && outerFlag4glbEl[ielGlb] == 1){
			lcl2glbEl[ielOuter] = ielGlb;
			glb2lclEl[ielGlb] = ielOuter;
			ielOuter++;
		}
	}
}
	
void countFaces(int nGlbElem, int rank, int *nFaces, int *neig, int *part, int *nGlbFaces, int *nLclFaces, int *nOuterCellFaces, int *nProxFaces, int *outerFlag4glbEl){
	
	int ielGlb1;

	int *seenElement = (int*)malloc(nGlbElem*sizeof(int));
	for(int ielGlb = 0; ielGlb < nGlbElem; ielGlb++){
		seenElement[ielGlb] = 0;
	}

	for(int ielGlb = 0; ielGlb < nGlbElem; ielGlb++){
		seenElement[ielGlb] = 1;
		for(int ifc = 0; ifc < nFaces[ielGlb]; ifc++){
			if(neig[ielGlb*6+ifc]==-1){
				ielGlb1=-1;
			}
			else{
				ielGlb1 = neig[ielGlb*6+ifc];
			}
			if(outerFlag4glbEl[ielGlb] == 0){
				if(ielGlb1 == -1 || seenElement[ielGlb1] == 0){
					(*nGlbFaces)++;
					if(ielGlb1 == -1 && part[ielGlb] == rank){
						(*nLclFaces)++;
					}else if(ielGlb1 > -1 && (part[ielGlb] == rank || part[ielGlb1] == rank)){
						(*nLclFaces)++;
					}
				}	
			}else{
				if(ielGlb1 == -1 || seenElement[ielGlb1] == 0){
					(*nGlbFaces)++;
					if(ielGlb1 == -1 && part[ielGlb] == rank){
						(*nLclFaces)++;
						(*nOuterCellFaces)++;	
					}else if(ielGlb1 > -1 && (part[ielGlb] == rank && part[ielGlb1] == rank)){
						(*nLclFaces)++;
						(*nOuterCellFaces)++;		
					}else if( ielGlb1 !=-1 && part[ielGlb] != rank && part[ielGlb1] == rank){
						nProxFaces[part[ielGlb]]++;
						(*nLclFaces)++;
						(*nOuterCellFaces)++;
					}else if(ielGlb1 !=-1 && part[ielGlb] == rank && part[ielGlb1] != rank){
						nProxFaces[part[ielGlb1]]++;
						(*nLclFaces)++;
						(*nOuterCellFaces)++;
					}
				}
			}	
		}	
	}
	free(seenElement);
}

void defFacePartitioning(int nGlbElem, int nLclFacesTot, int nProxTot, int nOuterCelllFaces, int rank, int *nFaces, int *neig, int *part, int *lcl2glbFc, int *glb2lclFc, int *outerFlag4glbEl){
	
	int ielGlb1;
	int nLcl  = 0;
	int nGlb  = 0;
	int idProxFace = nLclFacesTot - nProxTot;
	int idOuterCelllFace = nLclFacesTot - nOuterCelllFaces;

	int *seenElement = (int*)malloc(nGlbElem*sizeof(int));
	for(int ielGlb = 0; ielGlb < nGlbElem; ielGlb++){
		seenElement[ielGlb] = 0;
	}

	for(int ielGlb = 0; ielGlb < nGlbElem; ielGlb++){
		seenElement[ielGlb] = 1;
//		printf("iel is %d\n", ielGlb);
		for(int ifc = 0; ifc < nFaces[ielGlb]; ifc++){
			if(neig[ielGlb*6+ifc]==-1){
				ielGlb1=-1;
			}
			else{
				ielGlb1 = neig[ielGlb*6+ifc];
			}
//			printf("\tifc %d Nieg %d outerFlag %d\n", ifc, ielGlb1, outerFlag4glbEl[ielGlb]);
			if(ielGlb1 == -1 || seenElement[ielGlb1] == 0){
				if(outerFlag4glbEl[ielGlb] == 0){
					if(ielGlb1 == -1 && part[ielGlb] == rank){
						lcl2glbFc[nLcl] = nGlb;
						glb2lclFc[nGlb] = nLcl;
//						printf("\t\tBOUNDARY CELL lcl2glbFc[%d] =  %d, glb2lclFc[%d] = %d\n", nLcl, nGlb, nGlb, nLcl);
						nLcl++;
					}else if(ielGlb1 > -1 && (part[ielGlb] == rank && part[ielGlb1] == rank)){
						lcl2glbFc[nLcl] = nGlb;
						glb2lclFc[nGlb] = nLcl;
						//printf("\t\tFULL INNER lcl2glbFc[%d] =  %d, glb2lclFc[%d] = %d\n", nLcl, nGlb, nGlb, nLcl);
						nLcl++;
					}
				}
				if(outerFlag4glbEl[ielGlb] == 1){
					if(ielGlb1 == -1 && part[ielGlb] == rank){
						lcl2glbFc[idOuterCelllFace] = nGlb;
						glb2lclFc[nGlb] = idOuterCelllFace;
						//printf("\t\tOUTER BOUNDARY lcl2glbFc[%d] =  %d, glb2lclFc[%d] = %d\n", idOuterCelllFace, nGlb, nGlb, idOuterCelllFace);
						idOuterCelllFace++;
					}else if(ielGlb1 > -1 && (part[ielGlb] == rank && part[ielGlb1] == rank)){
						lcl2glbFc[idOuterCelllFace] = nGlb;
						glb2lclFc[nGlb] = idOuterCelllFace;
						//printf("\t\tOUTER NORMAL lcl2glbFc[%d] =  %d, glb2lclFc[%d] = %d\n", idOuterCelllFace, nGlb, nGlb, idOuterCelllFace);
						idOuterCelllFace++;
					}	
					else if(ielGlb1 > -1 && ( (part[ielGlb] == rank && part[ielGlb1] != rank) || (part[ielGlb] != rank && part[ielGlb1] == rank) ) ){
						lcl2glbFc[idProxFace] = nGlb;
						glb2lclFc[nGlb] = idProxFace;
						//printf("\t\tOUTER PROX lcl2glbFc[%d] =  %d, glb2lclFc[%d] = %d\n", idProxFace, nGlb, nGlb, idProxFace);
						idProxFace++;
					}
				}
				nGlb++;
			}
		}
	}
	free(seenElement);
}

void ghostCellsAllocation(int nGlbElem, int nLclElem, int rank, int size, int *nNeigRanks, int *lcl2glbEl, int *glb2lclEl, int *nFaces, int *neig, int *part, int **neigRanks){
	
	int ielGlb, ielGlb1, iel1;	
	int nGhostCells = 0;
	int iv =0;
	int bcFlag = 0; int proxFlag = 0;
	
	int *seenGhostElement = (int*)malloc(nGlbElem*sizeof(int));
	for(int iel = 0; iel < nGlbElem; iel++){
		seenGhostElement[iel] = 0;
	}

	int *seenRank = (int*)malloc(size*sizeof(int));
	int *neigRanks_temp = (int*)malloc(size*sizeof(int));
	for(int irank = 0; irank < size; irank++){
		seenRank[irank] = 0;
		neigRanks_temp[irank] = -1;
	}
	
	for(int iel = 0; iel < nLclElem; iel++){
		ielGlb = lcl2glbEl[iel];
		for(int ifc = 0; ifc < nFaces[ielGlb]; ifc++){
			bcFlag   = 0;
			proxFlag = 0;
			
			ielGlb1  = neig[ielGlb*6+ifc];
			bcFlag   = (ielGlb1 == -1);
			if(bcFlag == 0){
				proxFlag = (part[ielGlb1] != rank);
				if(proxFlag == 1){
					if(seenGhostElement[ielGlb1] == 0){
						iel1 = nLclElem + nGhostCells;
						glb2lclEl[ielGlb1] = iel1;
						lcl2glbEl[iel1]  = ielGlb1;
						seenGhostElement[ielGlb1] = 1;
						nGhostCells++;
					}
					if(!seenRank[part[ielGlb1]]){ 
						neigRanks_temp[iv] = part[ielGlb1];
						seenRank[part[ielGlb1]] = 1;
						iv++;
					}
				}
			}
		}
	}
	(*nNeigRanks) = iv;

	(*neigRanks) = (int*)malloc(iv*sizeof(int));
	for(int i = 0; i < iv; i++){
		(*neigRanks)[i] =  neigRanks_temp[i];
	}
	
	free(neigRanks_temp);
	free(seenGhostElement);
	free(seenRank);
}


void faceConnectivity(int nGlbElem, int nLclElem, int nLclFaces, int nOuterCellFaces, int nProxTot, int nGhostCells, int rank, int *lcl2glbEl, int *glb2lclEl, int *nFaces, int *neig, int *neigFc, int *part, int *fc2elFc, int *fc2el, int *elFc2fc, int *outerFlag4glbEl){

	int iel, ielGlb1, iel1, ifc1;	
	int bcFlag, proxFlag;
	int idOuter = nLclFaces - nOuterCellFaces;
	int idPrx   = nLclFaces - nProxTot;
	int idInner = 0;
	int countFlag = 0;

	
	for(int id = 0; id < nLclFaces; id++){
		fc2elFc[id] = -1;
		fc2el[id]  = -1;
	}
	for(int id = 0; id < nLclElem*6; id++){
		elFc2fc[id] = -1;
	}

	int *seenElement = (int*)malloc((nLclElem+nGhostCells)*sizeof(int));
	for(int iel = 0; iel < nLclElem+nGhostCells; iel++){
		seenElement[iel] = 0;
	}	
	
	for(int ielGlb = 0; ielGlb < nGlbElem; ielGlb++){
		if(glb2lclEl[ielGlb] > -1 ){
			iel = glb2lclEl[ielGlb];
			seenElement[iel] = 1;
			for(int ifc = 0; ifc < nFaces[ielGlb]; ifc++){
				bcFlag   = 0;
				proxFlag = 0;
				countFlag = 0;
				ifc1     = neigFc[ielGlb*6+ifc];
			
				ielGlb1  = neig[ielGlb*6+ifc];
				bcFlag   = (ielGlb1 == -1);
				if(bcFlag == 1){ 
					iel1 = -1;
					countFlag = (part[ielGlb] == rank);
				}
				else{
					iel1 = glb2lclEl[ielGlb1];
					countFlag = (part[ielGlb1] == rank || part[ielGlb] == rank);
					proxFlag = ( (part[ielGlb1] != rank && part[ielGlb] == rank) || (part[ielGlb1] == rank && part[ielGlb] != rank ) );
				}
				if(countFlag == 1){
				if(outerFlag4glbEl[ielGlb] == 0){
					if( iel1 == -1  || ( seenElement[iel1] == 0 && proxFlag != 1) ){
								
						elFc2fc[iel*6+ifc]   = idInner;
						if(iel1 != -1 && proxFlag != 1){elFc2fc[iel1*6+ifc1] = idInner;}
						
						fc2elFc[idInner*2 + 0] = ifc;
						fc2elFc[idInner*2 + 1] = ifc1;
					
						fc2el[idInner*2 + 0] = iel;
						fc2el[idInner*2 + 1] = iel1;

		
						idInner++;
					}
				}else{
					if(proxFlag == 1 && seenElement[iel1] == 0){
						elFc2fc[iel*6+ifc]   = idPrx;
						elFc2fc[iel1*6+ifc1] = idPrx;
				
						if(iel1 > iel){

							fc2elFc[idPrx*2 + 0] = ifc;
							fc2elFc[idPrx*2 + 1] = ifc1;
					
							fc2el[idPrx*2 + 0] = iel;
							fc2el[idPrx*2 + 1] = iel1;
						}else{
							fc2elFc[idPrx*2 + 0] = ifc1;
							fc2elFc[idPrx*2 + 1] = ifc;
					
							fc2el[idPrx*2 + 0] = iel1;
							fc2el[idPrx*2 + 1] = iel;
						}
	
						idPrx++;
					}else if(iel1 == -1  ||  seenElement[iel1] == 0){
						elFc2fc[iel*6+ifc]   = idOuter;
						if(iel1 != -1 ){elFc2fc[iel1*6+ifc1] = idOuter;}
				
						fc2elFc[idOuter*2 + 0] = ifc;
						fc2elFc[idOuter*2 + 1] = ifc1;
				
						fc2el[idOuter*2 + 0] = iel;
						fc2el[idOuter*2 + 1] = iel1;

						idOuter++;
					}
				}
				}
			}
		}
	}
	free(seenElement);
}


void localizeFaceValues(HOSTLOCALDATA *hLclData, int *lcl2glbEl, int *fc2el, int *fc2elFc, int *boundCond, double *fcCenter){
	
	int ielGlb, elFc;

	hLclData->boundCond = (int*)malloc(hLclData->nFcs*sizeof(int));
	hLclData->fcCenter  = (double*)malloc(hLclData->nFcs*3*sizeof(double));
	for(int ifc = 0; ifc < hLclData->nFcs; ifc++){
		ielGlb = lcl2glbEl[fc2el[ifc*2+0]];
		elFc   = fc2elFc[ifc*2+0];
		
		hLclData->boundCond[ifc] = boundCond[ielGlb*6 + elFc];

		for(int idr = 0; idr < 3; idr++){
			hLclData->fcCenter[ifc*3 + idr] = fcCenter[ielGlb*6*3 + elFc*3 + idr]; 
		}
	}
}


void countLclNodes(HOSTLOCALDATA *hLclData, int rank, int nGlbNodes, int max, int *nEl4GlbNd, int *idGlbEl4GlbNd,  int *part, int *glb2lclNd){

	int ielGlb;
	hLclData->nNodes = 0;
	for(int ignd = 0; ignd < nGlbNodes; ignd++){
		glb2lclNd[ignd] = -1;
		for(int iv = 0; iv < nEl4GlbNd[ignd]; iv++){
			ielGlb = idGlbEl4GlbNd[ignd*max + iv];
			if(part[ielGlb] == rank){
				glb2lclNd[ignd] = hLclData->nNodes; 
				hLclData->nNodes++;
				iv = 2e6;
			}
		}
	}
}


void localizeNodeValues(HOSTLOCALDATA *hLclData, int nGlbNodes, int max, int *glb2lclEl, int *nEl4GlbNd, int *idGlbEl4GlbNd, int *glb2lclNd, double *glbNdCord){
				
	int ilnd;

	hLclData->lcl2glbNd = (int*)malloc(hLclData->nNodes*sizeof(int));
	hLclData->nEl4nd    = (int*)malloc(hLclData->nNodes*sizeof(int));
	hLclData->idEl4nd   = (int*)malloc(hLclData->nNodes*max*sizeof(int));
	hLclData->glbNdCord = (double*)malloc(hLclData->nNodes*3*sizeof(double));
	for(int ignd = 0; ignd < nGlbNodes; ignd++){
		ilnd = glb2lclNd[ignd];  
		if(ilnd!=-1){
			hLclData->lcl2glbNd[ilnd] = ignd;
			hLclData->nEl4nd[ilnd] = nEl4GlbNd[ignd];
			for(int iv = 0; iv < nEl4GlbNd[ignd]; iv++){
				hLclData->idEl4nd[ilnd*max + iv] = glb2lclEl[idGlbEl4GlbNd[ignd*max + iv]];
				if(glb2lclEl[idGlbEl4GlbNd[ignd*max + iv]] > hLclData->nElem-1){ // Do not include ghost cells
					hLclData->idEl4nd[ilnd*max + iv] = -1;
				}
			}
			for(int idr = 0; idr < 3; idr++){
				hLclData->glbNdCord[ilnd*3+idr] = glbNdCord[ignd*3+idr];
			}
		}
	}
}


void localizeCellCenter(HOSTLOCALDATA *hLclData, int *lcl2glbEl, int *glbElNd2GlbNd, int *glb2lclNd, double *cellCenter, int *elNd2lclNd, int *nLclNodes, int *nLclNodes4glbEl){

	int ielGlb;
	hLclData->cellCenter = (double*)malloc(hLclData->nElem*3*sizeof(double));
	for(int iel = 0; iel < hLclData->nElem; iel++){
		ielGlb = lcl2glbEl[iel];
		nLclNodes[iel] = nLclNodes4glbEl[ielGlb];
		for(int idr = 0; idr < 3; idr++){
			hLclData->cellCenter[iel*3 + idr] = cellCenter[ielGlb*3 + idr];
		}
		for(int ind = 0; ind < 8; ind++){
			elNd2lclNd[iel*8+ind] = glb2lclNd[glbElNd2GlbNd[ielGlb*8 + ind]];
		}
	}
}
	
void communicateVolumes(HOSTLOCALDATA *hLclData, int rank, int size, int nNeigRanks, int max, int *neigRanks, int *nProxFaces, int *neigRank4fc, int *lclProx4fc, int * lclFc2idRcv, int *fc2el, double *volume){ 
	
	int neigRank, iProx;
	double *sendBuff = (double*)malloc(size*max*sizeof(double));
	double *recvBuff = (double*)malloc(size*max*sizeof(double));
	
	for(int ifc = hLclData->nFcs - hLclData->nProxTot; ifc < hLclData->nFcs; ifc++){
		neigRank = neigRank4fc[ifc];
		iProx    = lclProx4fc[ifc];
		sendBuff[neigRank*max + iProx] = volume[fc2el[ifc*2+0]]; 
	}
	
	MPI_Request reqs[2*nNeigRanks];
	int r=0;
	for(int iv = 0; iv < nNeigRanks; iv++){
		neigRank = neigRanks[iv];
		MPI_Isend(&sendBuff[neigRank*max],nProxFaces[neigRank],MPI_DOUBLE,neigRank,rank,MPI_COMM_WORLD,&reqs[r++]);
	}	
	for(int ipt = 0; ipt < nNeigRanks; ipt++){
		neigRank = neigRanks[ipt];
		MPI_Irecv(&recvBuff[neigRank*max],nProxFaces[neigRank],MPI_DOUBLE,neigRank,neigRank,MPI_COMM_WORLD,&reqs[r++]);
	}
	MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);

	for(int ifc = hLclData->nFcs - hLclData->nProxTot; ifc < hLclData->nFcs; ifc++){
		neigRank = neigRank4fc[ifc];
		iProx    = lclFc2idRcv[ifc];
		volume[fc2el[ifc*2+1]] = recvBuff[neigRank*max + iProx];
	}
	free(sendBuff);
	free(recvBuff);
}

        
void buildLocalMesh(HOSTLOCALDATA *hLclData, int nGlbElem, int nGlbNodes, int max, int *nLclNodes4glbEl, int *neig, int *neigFc, int *nFaces, int *part, int *boundCond, int *nEl4GlbNd, int *idGlbEl4GlbNd, int *glbElNd2GlbNd, double *cellCenter, double *fcCenter, double *glbNdCord){
	
	int rank, size;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	int nLclElem    = countLocalElem(nGlbElem,rank,part);
	int nOuter;
	int *outerFlag4glbEl = (int*)malloc(nGlbElem*sizeof(int));
	int nGhostCells = countGhostCells(nGlbElem, nLclElem, rank, &nOuter, nFaces, part, neig, outerFlag4glbEl); 
	MPI_Allreduce(MPI_IN_PLACE, outerFlag4glbEl, nGlbElem, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	int *glb2lclEl = (int*)malloc(nGlbElem*sizeof(int));
	int *lcl2glbEl = (int*)malloc((nLclElem+nGhostCells)*sizeof(int));
	defElemPartitioning(nGlbElem, nLclElem, nOuter, rank, part, lcl2glbEl, glb2lclEl, outerFlag4glbEl);
	
	int nNeigRanks = 0;
	int *neigRanks;
	ghostCellsAllocation(nGlbElem, nLclElem, rank, size, &nNeigRanks, lcl2glbEl, glb2lclEl, nFaces, neig, part, &neigRanks);

	int nGlbFaces = 0; int nLclFaces = 0; int nOuterCellFaces =0;
	int *nProxFaces = (int*)malloc(size*sizeof(int));
	for(int irank = 0; irank < size; irank++){
		nProxFaces[irank]=0;
	}
	countFaces(nGlbElem, rank, nFaces, neig, part, &nGlbFaces, &nLclFaces, &nOuterCellFaces, nProxFaces, outerFlag4glbEl);
	
	if(rank==0){printf("Faces per rank: ");}
	for(int irank = 0; irank < size; irank++){
		MPI_Barrier(MPI_COMM_WORLD);
		if(irank==rank){printf("%d ", nLclFaces);fflush(stdout);}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){printf("\n");}

	int *glb2lclFc = (int*)malloc(nGlbFaces*sizeof(int));
	for(int iGlbFc = 0; iGlbFc < nGlbFaces; iGlbFc++){
		glb2lclFc[iGlbFc] = -1;
	}
	int *lcl2glbFc = (int*)malloc(nLclFaces*sizeof(int));
	int nProxTot = 0;
	for(int irank = 0; irank < size; irank++){
		nProxTot += nProxFaces[irank];
	}
	defFacePartitioning(nGlbElem, nLclFaces, nProxTot, nOuterCellFaces, rank, nFaces, neig, part, lcl2glbFc, glb2lclFc, outerFlag4glbEl);
	
	int *fc2elFc = (int*)malloc(nLclFaces*2*sizeof(int));
	int *fc2el   = (int*)malloc(nLclFaces*2*sizeof(int));
	int *elFc2fc = (int*)malloc((nLclElem+nGhostCells)*6*sizeof(int));
	faceConnectivity(nGlbElem, nLclElem, nLclFaces, nOuterCellFaces, nProxTot, nGhostCells, rank, lcl2glbEl, glb2lclEl, nFaces, neig, neigFc, part, fc2elFc, fc2el, elFc2fc, outerFlag4glbEl);
	free(neigFc); free(outerFlag4glbEl);
		
	int nProxFacesMax = 0;
	int irnk;
	for(int ipr = 0; ipr < nNeigRanks; ipr++){
		irnk = neigRanks[ipr];
		if(nProxFaces[irnk]>nProxFacesMax){
			nProxFacesMax = nProxFaces[irnk];
		}
	}

	hLclData->neigRank4fc = (int*)malloc(nLclFaces*sizeof(int));
	hLclData->lclProx4fc  = (int*)malloc(nLclFaces*sizeof(int));
	hLclData->lclFc2idRcv = (int*)malloc(nLclFaces*sizeof(int));
	for(int i = 0; i < nLclFaces; i++){
		hLclData->neigRank4fc[i] = -1;
		hLclData->lclProx4fc[i] = -1;
		hLclData->lclFc2idRcv[i] = -1;
	}
	partitionConnectivity(nLclElem, rank, size, nNeigRanks, nProxFacesMax, nFaces, lcl2glbEl, lcl2glbFc, glb2lclFc, neig, part, elFc2fc, neigRanks, nProxFaces, hLclData->neigRank4fc, hLclData->lclProx4fc, hLclData->lclFc2idRcv);
	free(neig); free(nFaces); free(glb2lclFc); free(lcl2glbFc); free(elFc2fc);
	
	hLclData->nElem      = nLclElem;
	hLclData->nGhostElem = nGhostCells;
	hLclData->nOuter     = nOuter;
	hLclData->nProxTot   = nProxTot;
	hLclData->nFcs       = nLclFaces;
	localizeFaceValues(hLclData, lcl2glbEl, fc2el, fc2elFc, boundCond, fcCenter);
	free(boundCond); free(fcCenter);

	int *glb2lclNd = (int*)malloc(nGlbNodes*sizeof(int));
	countLclNodes(hLclData, rank, nGlbNodes, max, nEl4GlbNd, idGlbEl4GlbNd, part, glb2lclNd);
	free(part);
	
	localizeNodeValues(hLclData, nGlbNodes, max, glb2lclEl, nEl4GlbNd, idGlbEl4GlbNd, glb2lclNd, glbNdCord);
	free(nEl4GlbNd); free(idGlbEl4GlbNd); free(glbNdCord); free(glb2lclEl);

	int *elNd2lclNd = (int*)malloc(hLclData->nElem*8*sizeof(int));
	int *nLclNodes  = (int*)malloc(hLclData->nElem*sizeof(int));
	localizeCellCenter(hLclData, lcl2glbEl, glbElNd2GlbNd, glb2lclNd, cellCenter, elNd2lclNd, nLclNodes, nLclNodes4glbEl);
	hLclData->maxEl4nd = max;
	free(cellCenter); free(glb2lclNd); free(glbElNd2GlbNd); free(nLclNodes4glbEl); 

	meshProperties(hLclData, nLclNodes, fc2el, fc2elFc, elNd2lclNd);
	free(fc2elFc);

	communicateVolumes(hLclData, rank, size, nNeigRanks, nProxFacesMax, neigRanks, nProxFaces, hLclData->neigRank4fc, hLclData->lclProx4fc, hLclData->lclFc2idRcv, fc2el, hLclData->volume);

	hLclData->nNeigRanks = nNeigRanks;
	hLclData->nProxFacesMax = nProxFacesMax;
	hLclData->nProxFaces = (int*)malloc(size*sizeof(int));
	for(int i = 0; i < size; i++){
		hLclData->nProxFaces[i] = nProxFaces[i];
	}
	free(nProxFaces);
	hLclData->neigRanks = (int*)malloc(nNeigRanks*sizeof(int));
	for(int i = 0; i < nNeigRanks; i++){
		hLclData->neigRanks[i] = neigRanks[i];
	}
	free(neigRanks);
	hLclData->fc2el = (int*)malloc(nLclFaces*2*sizeof(int));
	for(int i = 0; i < nLclFaces*2; i++){
		hLclData->fc2el[i] = fc2el[i];
	}
	free(fc2el);
	hLclData->lcl2glbEl = (int*)malloc(nLclElem*sizeof(int));
	for(int i = 0; i < nLclElem; i++){
		hLclData->lcl2glbEl[i] = lcl2glbEl[i];
	}
	free(lcl2glbEl);
	hLclData->nElNds = (int*)malloc(nLclElem*sizeof(int));
	for(int i = 0; i < nLclElem; i++){
		hLclData->nElNds[i] = nLclNodes[i];
	}
	free(nLclNodes);
	hLclData->elNd2lclNd = (int*)malloc(nLclElem*8*sizeof(int));
	for(int i = 0; i < nLclElem*8; i++){
		hLclData->elNd2lclNd[i] = elNd2lclNd[i];
	}
	free(elNd2lclNd);	
	
}

//	freeVars();

//	if(rank==3){
//		printf("lcl2glbNd[3] %d nEl4nd[3] %d idEl4nd[3][0] %d idEl4nd[3][1] %d x %f y %f\n", hLclData->lcl2glbNd[0], hLclData->nEl4nd[0], hLclData->idEl4nd[0*max + 0], hLclData->idEl4nd[0*max+1], hLclData->glbNdCord[0*3 + 0],  hLclData->glbNdCord[0*3 + 1]);
//	}

//	if(rank==2){
//		printf("ifc 12 xc yc %lf %lf\n", hLclData->fcCenter[12*3+0], hLclData->fcCenter[12*3+1]);
//		printf("ifc 21 xc yc %lf %lf\n", hLclData->fcCenter[21*3+0], hLclData->fcCenter[21*3+1]);
//		printf("ifc 38 xc yc %lf %lf\n", hLclData->fcCenter[38*3+0], hLclData->fcCenter[38*3+1]);
//		printf("ifc 48 xc yc %lf %lf\n", hLclData->fcCenter[48*3+0], hLclData->fcCenter[48*3+1]);
//		printf("ifc 18 xc yc %lf %lf\n", hLclData->fcCenter[18*3+0], hLclData->fcCenter[18*3+1]);
//		printf("ifc 12 bc is %d\n", hLclData->boundCond[12]);
//		printf("ifc 18 bc is %d\n", hLclData->boundCond[18]);
//		printf("ifc 6  bc is %d\n", hLclData->boundCond[6]);
//		printf("ifc 23 bc is %d\n", hLclData->boundCond[23]);
//		printf("max %d ind 26 %d\n",max, h_idGlbEl4GlbNd[26*max + 0]);
	//	printf("nNodes %d\n", hLclData->nNodes);
	//	printf("lcl nd 6 x y %f %f\n", hLclData->glbNdCord[6*3 + 0], hLclData->glbNdCord[6*3 + 1]); 
	//	printf("lcl nd 11 elements %d %d %d %d\n", hLclData->idEl4nd[11*max + 0], hLclData->idEl4nd[11*max + 1],hLclData->idEl4nd[11*max + 2],hLclData->idEl4nd[11*max + 3]);
	//	printf("cell 4 %lf %lf\n", hLclData->cellCenter[4*3+0], hLclData->cellCenter[4*3+1]);
	//	printf("cell 8 %lf %lf\n", hLclData->cellCenter[8*3+0], hLclData->cellCenter[8*3+1]);
//		printf("Neig Ranks %d %d\n", hLclData->neigRanks[0], hLclData->neigRanks[1]);
	//	printf("nElem %d nGhostCells %d nGlbFaces %d nLclFaces %d nProxTot %d nOuterCellsFaces %d\n", nLclElem, nGhostCells, nGlbFaces, nLclFaces, nProxTot, nOuterCellFaces);
//		printf("Prox faces with rank 0 %d rank 1 %d rank 3 %d\n", nProxFaces[0], nProxFaces[1],nProxFaces[3]);
// 		printf("lcl2glbEl[4] %d lcl2glbEl[10] %d glb2lcl[6] %d glb2lcl[4] %d\n", lcl2glbEl[4], lcl2glbEl[10], glb2lclEl[6], glb2lclEl[4]); 
//		printf("lcl2glbFc[49] %d lcl2glbFc[53] %d glb2lclFc[12] %d glb2lclFc[86] %d glb2lclFc[64] %d\n", lcl2glbFc[49], lcl2glbFc[53], glb2lclFc[12], glb2lclFc[86], glb2lclFc[64]);	
	//	printf("lcl2glbFc[48] %d lcl2glbFc[49] %d lcl2glbFc[50] %d lcl2glbFc[51] %d lcl2glbFc[52] %d lcl2glbFc[53] %d lcl2glbFc[54] %d\n",lcl2glbFc[48],lcl2glbFc[49],lcl2glbFc[50],lcl2glbFc[51],lcl2glbFc[52],lcl2glbFc[53],lcl2glbFc[54]);
//		printf("fc2el[17] %d %d fc2el[36] %d %d fc2el[22] %d %d fc2el[53] %d %d\n", fc2el[17*2+0], fc2el[17*2+1], fc2el[36*2+0], fc2el[36*2+1], fc2el[22*2+0], fc2el[22*2+1],fc2el[53*2+0], fc2el[53*2+1]); 	
//		printf("fc2elFc[17] %d %d fc2elFc[36] %d %d fc2elFc[22] %d %d fc2elFc[53] %d %d\n", fc2elFc[17*2+0], fc2elFc[17*2+1], fc2elFc[36*2+0], fc2elFc[36*2+1], fc2elFc[22*2+0], fc2elFc[22*2+1],fc2elFc[53*2+0], fc2elFc[53*2+1]); 	
//		printf("elFc2fc[3][1] %d elFc2fc[10][2] %d elFc2fc[4][3] %d elFc2fc[13][3] %d\n", elFc2fc[3*6 + 1], elFc2fc[10*6 + 2],elFc2fc[4*6 + 3],elFc2fc[13*6 + 3]);

//	}
	
//	if(rank==3){
//		printf("nElem %d nGhostCells %d nGlbFaces %d nLclFaces %d nProxTot %d nOuterCellsFaces %d\n", nLclElem, nGhostCells, nGlbFaces, nLclFaces, nProxTot, nOuterCellFaces);
//		printf("Prox faces with rank 0 %d rank 1 %d rank 3 %d\n", nProxFaces[0], nProxFaces[1],nProxFaces[3]);
//		printf("lcl2glbFc[48] %d lcl2glbFc[49] %d lcl2glbFc[50] %d lcl2glbFc[51] %d lcl2glbFc[52] %d lcl2glbFc[53] %d lcl2glbFc[54] %d\n",lcl2glbFc[48],lcl2glbFc[49],lcl2glbFc[50],lcl2glbFc[51],lcl2glbFc[52],lcl2glbFc[53],lcl2glbFc[54]);
	//	printf("glb2lclEl[2] %d glb2lclEl[8] %d glb2lclEl[14] %d glb2lclEl[20] %d glb2lclEl[27] %d glb2lclEl[28] %d glb2lclEl[29] %d\n",glb2lclEl[2],glb2lclEl[8],glb2lclEl[14],glb2lclEl[20],glb2lclEl[27],glb2lclEl[28],glb2lclEl[29]);
//	}
