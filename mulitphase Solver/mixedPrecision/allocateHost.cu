#include "strdata.h"

void allocateHost(int nGlbElem, int nGlbNodes, int **h_neig, int **h_neigFc,int **h_nEl4GlbNd, int **h_part){

	(*h_neig)      = (int*)malloc(nGlbElem*6*sizeof(int));
	(*h_neigFc)    = (int*)malloc(nGlbElem*6*sizeof(int));
	(*h_part)      = (int*)malloc(nGlbElem*sizeof(int));
	(*h_nEl4GlbNd) = (int*)malloc(nGlbNodes*sizeof(int));

}
