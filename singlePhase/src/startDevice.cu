#include "strdata.h"

void startDevice(int nGlbElem, int nGlbNodes, int rank, int size, HOSTLOCALDATA *hLclData, MESH *mesh, FIELD *field, COMM *comm, EXPORT *exp){

	int nFcs          = hLclData->nFcs;
	int nElem         = hLclData->nElem;
	int nGhostElem    = hLclData->nGhostElem;
	int nNeigRanks    = hLclData->nNeigRanks;
	int nProxFacesMax = hLclData->nProxFacesMax;
	int nNodes        = hLclData->nNodes;
	int maxEl4nd      = hLclData->maxEl4nd;

	size_t offset = 0;
	unsigned char *d_buffer;

	// Mesh Struct
	size_t bytes_fc2el      = nFcs  	       * 2 * sizeof(int);
	size_t bytes_boundCond  = nFcs                 *     sizeof(int);
	size_t bytes_cellCenter = nElem                * 3 * sizeof(float);
	size_t bytes_fcCenter   = nFcs  	       * 3 * sizeof(float);
	size_t bytes_n          = nFcs  	       * 3 * sizeof(float);
	size_t bytes_nt1        = nFcs  	       * 3 * sizeof(float);
	size_t bytes_nt2        = nFcs  	       * 3 * sizeof(float);
	size_t bytes_area       = nFcs  	       *     sizeof(float);
	size_t bytes_volume     = (nElem + nGhostElem) *     sizeof(float);
	
	size_t meshTotalBytes = bytes_fc2el    + bytes_boundCond + bytes_cellCenter +
				bytes_fcCenter + bytes_n         + bytes_nt1 +	
				bytes_nt2      + bytes_area      + bytes_volume;
	//Field Struct
	size_t bytes_consVec0     = nElem                * 5 *     sizeof(float);
	size_t bytes_consVec      = (nElem + nGhostElem) * 5 *     sizeof(float);
	size_t bytes_primVecF     = nFcs   		 * 2 * 5 * sizeof(float);
	size_t bytes_RHS          = (nElem + nGhostElem) * 5 *     sizeof(float);
	size_t bytes_grad         = (nElem + nGhostElem) * 5 * 3 * sizeof(float);
	size_t bytes_wMax         = (nElem + nGhostElem) * 5 *     sizeof(float);
	size_t bytes_wMin         = (nElem + nGhostElem) * 5 *     sizeof(float);
	size_t bytes_extensionVec = nFcs  		 * 2 * 5 * sizeof(float);
	size_t bytes_theta        = (nElem + nGhostElem) * 5 *     sizeof(float);

	size_t fieldTotalBytes = bytes_consVec0 + bytes_consVec + bytes_primVecF     + bytes_RHS + bytes_grad +
				 bytes_wMax     + bytes_wMin    + bytes_extensionVec + bytes_theta;

	// Communication Struct -> (They can be reduce to nNeigRanks instead of size !)
	size_t bytes_neigRank4fc        = nFcs *                         sizeof(int);  
	size_t bytes_lclProx4fc         = nFcs *                         sizeof(int);  
	size_t bytes_lclFc2idRcv        = nFcs *                         sizeof(int);  
	size_t bytes_sendbuff           = size * nProxFacesMax * 5 *     sizeof(float);
	size_t bytes_recvbuff           = size * nProxFacesMax * 5 *     sizeof(float);
	size_t bytes_sendBuffGrad       = size * nProxFacesMax * 5 * 3 * sizeof(float);
	size_t bytes_recvBuffGrad       = size * nProxFacesMax * 5 * 3 * sizeof(float);

	size_t commTotalBytes = bytes_neigRank4fc + bytes_lclProx4fc + bytes_lclFc2idRcv  +   
			        bytes_sendbuff    + bytes_recvbuff   + bytes_sendBuffGrad +
				bytes_recvBuffGrad;

	// Export struct
	size_t bytes_lcl2glbEl  = nElem  *            sizeof(int);
	size_t bytes_lcl2glbNd  = nNodes *            sizeof(int);
	size_t bytes_nElNds     = nElem  *            sizeof(int);
	size_t bytes_elNd2lclNd = nElem  * 8        * sizeof(int);
	size_t bytes_nEl4nd     = nNodes *            sizeof(int);
	size_t bytes_idEl4nd    = nNodes * maxEl4nd * sizeof(int);
	size_t bytes_glbNdCord  = nNodes * 3        * sizeof(float); 

	size_t expTotalBytes = bytes_lcl2glbEl + bytes_lcl2glbNd + bytes_nElNds   + bytes_elNd2lclNd +
			       bytes_nEl4nd    + bytes_idEl4nd   + bytes_glbNdCord;

	size_t totalBytes    = meshTotalBytes + fieldTotalBytes + commTotalBytes + expTotalBytes;

	cudaMalloc(&d_buffer, totalBytes);

	mesh->fc2el      = (int*)(d_buffer + offset);
	offset += bytes_fc2el;

	mesh->boundCond  = (int*)(d_buffer + offset);
	offset += bytes_boundCond;

	mesh->cellCenter = (float*)(d_buffer + offset);
	offset += bytes_cellCenter;

	mesh->fcCenter   = (float*)(d_buffer + offset);
	offset += bytes_fcCenter;

	mesh->n          = (float*)(d_buffer + offset);
	offset += bytes_n;

	mesh->nt1        = (float*)(d_buffer + offset);
	offset += bytes_nt1;

	mesh->nt2        = (float*)(d_buffer + offset);
	offset += bytes_nt2;

	mesh->area       = (float*)(d_buffer + offset);
	offset += bytes_area;

	mesh->volume     = (float*)(d_buffer + offset);
	offset += bytes_volume;
	
	field->consVec0 = (float*)(d_buffer + offset);
	offset += bytes_consVec0;

	field->consVec = (float*)(d_buffer + offset);
	offset += bytes_consVec;

	field->primVecF = (float*)(d_buffer + offset);
	offset += bytes_primVecF;

	field->RHS = (float*)(d_buffer + offset);
	offset += bytes_RHS;

	field->grad = (float*)(d_buffer + offset);
	offset += bytes_grad;

	field->wMax = (float*)(d_buffer + offset);
	offset += bytes_wMax;

	field->wMin = (float*)(d_buffer + offset);
	offset += bytes_wMin;

	field->extensionVec = (float*)(d_buffer + offset);
	offset += bytes_extensionVec;

	field->theta = (float*)(d_buffer + offset);
	offset += bytes_theta;

	comm->neigRank4fc = (int*)(d_buffer + offset);
	offset += bytes_neigRank4fc;

	comm->lclProx4fc = (int*)(d_buffer + offset);
	offset += bytes_lclProx4fc;

	comm->lclFc2idRcv = (int*)(d_buffer + offset);
	offset += bytes_lclFc2idRcv;

	comm->sendbuff = (float*)(d_buffer + offset);
	offset += bytes_sendbuff;

	comm->recvbuff = (float*)(d_buffer + offset);
	offset += bytes_recvbuff;

	comm->sendBuffGrad = (float*)(d_buffer + offset);
	offset += bytes_sendBuffGrad;
	
	comm->recvBuffGrad = (float*)(d_buffer + offset);
	offset += bytes_recvBuffGrad;

	exp->lcl2glbEl = (int*)(d_buffer + offset);
	offset += bytes_lcl2glbEl;

	exp->lcl2glbNd = (int*)(d_buffer + offset);
	offset += bytes_lcl2glbNd;

	exp->nElNds = (int*)(d_buffer + offset);
	offset += bytes_nElNds;

	exp->elNd2lclNd = (int*)(d_buffer + offset);
	offset += bytes_elNd2lclNd;

	exp->nEl4nd = (int*)(d_buffer + offset);
	offset += bytes_nEl4nd;

	exp->idEl4nd = (int*)(d_buffer + offset);
	offset += bytes_idEl4nd;

	exp->glbNdCord = (float*)(d_buffer + offset);
	offset += bytes_glbNdCord;

	cudaMemcpy(mesh->fc2el,       hLclData->fc2el,       nFcs * 2 * sizeof(int),    cudaMemcpyHostToDevice);
	free(hLclData->fc2el);

	cudaMemcpy(mesh->boundCond,   hLclData->boundCond,   nFcs * sizeof(int),        cudaMemcpyHostToDevice);
	free(hLclData->boundCond);

	cudaMemcpy(mesh->cellCenter,  hLclData->cellCenter,  nElem * 3 * sizeof(float), cudaMemcpyHostToDevice);
	free(hLclData->cellCenter);

	cudaMemcpy(mesh->fcCenter,    hLclData->fcCenter,    nFcs  * 3 * sizeof(float), cudaMemcpyHostToDevice);
	free(hLclData->fcCenter);

	cudaMemcpy(mesh->n,           hLclData->n,           nFcs  * 3 * sizeof(float), cudaMemcpyHostToDevice);
	free(hLclData->n);

	cudaMemcpy(mesh->nt1,         hLclData->nt1,         nFcs  * 3 * sizeof(float), cudaMemcpyHostToDevice);
	free(hLclData->nt1);

	cudaMemcpy(mesh->nt2,         hLclData->nt2,         nFcs  * 3 * sizeof(float), cudaMemcpyHostToDevice);
	free(hLclData->nt2);

	cudaMemcpy(mesh->area,        hLclData->area,        nFcs  * sizeof(float),     cudaMemcpyHostToDevice);
	free(hLclData->area);

	cudaMemcpy(mesh->volume,      hLclData->volume,      (nElem + nGhostElem) * sizeof(float),     cudaMemcpyHostToDevice);
	free(hLclData->volume);

	cudaMemcpy(comm->neigRank4fc, hLclData->neigRank4fc,
	           nFcs * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->neigRank4fc);

	cudaMemcpy(comm->lclProx4fc, hLclData->lclProx4fc,
	           nFcs * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lclProx4fc);
	
	cudaMemcpy(comm->lclFc2idRcv, hLclData->lclFc2idRcv,
	           nFcs * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lclFc2idRcv);

	cudaMemcpy(exp->lcl2glbEl, hLclData->lcl2glbEl,
	           nElem * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lcl2glbEl);

	cudaMemcpy(exp->lcl2glbNd, hLclData->lcl2glbNd,
	           nNodes * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->lcl2glbNd);

	cudaMemcpy(exp->nElNds, hLclData->nElNds,
		  nElem * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->nElNds);

	cudaMemcpy(exp->elNd2lclNd, hLclData->elNd2lclNd,
		   nElem * 8 * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->elNd2lclNd);

	cudaMemcpy(exp->nEl4nd, hLclData->nEl4nd,
	           nNodes * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->nEl4nd);

	cudaMemcpy(exp->idEl4nd, hLclData->idEl4nd,
	           nNodes * maxEl4nd * sizeof(int), cudaMemcpyHostToDevice);
	free(hLclData->idEl4nd);

	cudaMemcpy(exp->glbNdCord, hLclData->glbNdCord,
	           nNodes * 3 * sizeof(float), cudaMemcpyHostToDevice);
	free(hLclData->glbNdCord);

	mesh->nElem         = nElem;
	mesh->nFcs          = nFcs;
	mesh->nOuter        = hLclData->nOuter;
	mesh->nProxTot      = hLclData->nProxTot;
	mesh->nGhostElem    = hLclData->nGhostElem;
	comm->nNeigRanks    = nNeigRanks;
	comm->nProxFacesMax = nProxFacesMax;
	comm->size          = size;
	comm->rank          = rank;
	exp->nGlbElem       = nGlbElem;
	exp->nGlbNodes      = nGlbNodes;
	exp->nNodes         = nNodes;
	exp->nElem          = nElem;
	exp->maxEl4nd       = maxEl4nd;
	
	comm->neigRanks  = (int*)malloc(nNeigRanks*sizeof(int));
	for(int i = 0; i < comm->nNeigRanks; i++){
		comm->neigRanks[i] = hLclData->neigRanks[i];
	}
	free(hLclData->neigRanks);
	
	comm->nProxFaces = (int*)malloc(nNeigRanks*sizeof(int));
	for(int i = 0; i < size; i++){
		comm->nProxFaces[i] = hLclData->nProxFaces[i];
	}
	free(hLclData->nProxFaces);

	comm->recvRequests = (MPI_Request**)malloc(3 * sizeof(MPI_Request *));
	for(int i = 0; i < 3; i++){
		comm->recvRequests[i] = (MPI_Request*)malloc(nNeigRanks * sizeof(MPI_Request));
	}
	
	free(hLclData);	 
	
}
	
