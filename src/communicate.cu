#include "strdata.h"

void sendArray(COMM * comm, int iField, int size, float *sendArr){
	int neigRank;
	for(int iv = 0; iv < comm->nNeigRanks; iv++){
		neigRank = comm->neigRanks[iv];
		MPI_Request request;
		MPI_Isend(&sendArr[neigRank*comm->nProxFacesMax*size],comm->nProxFaces[neigRank]*size,MPI_FLOAT,neigRank,(comm->rank+(2*comm->size*iField)),MPI_COMM_WORLD,&request);
		MPI_Request_free(&request);
	}

}

void recvArray(COMM * comm, int iField, int size, float *arr){

	int neigRank;
	for (int ipt = 0; ipt < comm->nNeigRanks; ipt++){
		neigRank = comm->neigRanks[ipt];
		MPI_Irecv(&arr[neigRank*comm->nProxFacesMax*size],comm->nProxFaces[neigRank]*size,MPI_FLOAT,neigRank,(neigRank+(2*comm->size*iField)),MPI_COMM_WORLD,&comm->recvRequests[iField][ipt]);
	}
}

