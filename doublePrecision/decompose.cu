#include "strdata.h"
	
void decompose(int nGlbElem, int *h_nFaces, int *h_neig, int *h_part){
	
	int *XADJ, *ADJY;
	int *vwgt;
	int iedge,ineig;
	int nedge, ncon;
	int objval,size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	XADJ = (int*)malloc((nGlbElem+1)*sizeof(int));
	vwgt = (int*)malloc((nGlbElem)  *sizeof(int));
		
	nedge = 0;
	for (int iel = 0; iel < nGlbElem; iel++){
		for (int ifc = 0; ifc < h_nFaces[iel]; ifc++){
			if (h_neig[iel*6 + ifc] >= 0) {
				nedge++; 
			}
		}
	}
	
	ADJY = (int*)malloc(nedge*sizeof(int));

	iedge=0;
	ineig=0;
	for (int iel = 0; iel < nGlbElem; iel++){
		vwgt[iel] = h_nFaces[iel];
		if (iel == 0) {
			XADJ[iel]=0;
		}
		if (iel != 0) {	
			XADJ[iel]=XADJ[iel-1]+ineig;
		}
		ineig=0;
		for (int ifc = 0; ifc < h_nFaces[iel]; ifc++){
			if (h_neig[iel*6 + ifc] >= 0) {
				ADJY[iedge] = h_neig[iel*6 + ifc];
				iedge++;
				ineig++;
			}
		}
	}
	XADJ[nGlbElem]=XADJ[nGlbElem-1]+ineig;
	
	ncon = 1;
	if(size > 1){
		METIS_PartGraphRecursive(&nGlbElem,&ncon,XADJ,ADJY,vwgt,NULL,NULL,&size,NULL,NULL,NULL,&objval,h_part); 
	}else{
		for(int iel = 0; iel < nGlbElem; iel++){
			h_part[iel] = 0;
		}
	}


//	if(rank==2){
//		for(int iel = 0; iel < nGlbElem; iel++){
//			printf("iel %d part %d\n", iel, h_part[iel]);
//		}
//	}
	
	free(vwgt);
	free(XADJ);
	free(ADJY);
}
