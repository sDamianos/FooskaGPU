#include "strdata.h"

void mainStart(MESH *mesh, EXPORT *exp, COMM *comm, FIELD *field){


	int nGlbNodes, nGlbElem, maxEl4GlbNd;
	
	//Host Allocations
	int   *h_nFaces, *h_nLclNodes, *h_elNd2GlbNd, *h_boundCond;
	int   *h_neig,   *h_neigFc, *h_part;
	int   *h_idGlbEl4GlbNd, *h_nEl4GlbNd;
	float *h_cellCenter, *h_fcCenter, *h_glbNdCrdnts;

	int rank,size;
	HOSTLOCALDATA *hLclData   = (HOSTLOCALDATA*) malloc(sizeof(HOSTLOCALDATA));
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	translator_ptw(input.scale,&nGlbElem,&nGlbNodes,&h_nFaces,&h_nLclNodes,&h_elNd2GlbNd,&h_boundCond,&h_cellCenter,&h_fcCenter,&h_glbNdCrdnts);
	if(rank == 0){printf("Neutral file...........: OK \n");}
	
	allocateHost(nGlbElem,nGlbNodes,&h_neig,&h_neigFc,&h_nEl4GlbNd,&h_part);
	
	meshConnectivity(nGlbElem,nGlbNodes,h_nFaces,h_nLclNodes,h_elNd2GlbNd,h_neig,h_neigFc,h_boundCond,h_nEl4GlbNd,&h_idGlbEl4GlbNd,&maxEl4GlbNd,h_fcCenter); 
	if(rank == 0){printf("Mesh Connectivity......: OK \n");}

	decompose(nGlbElem,h_nFaces,h_neig,h_part);
	if(rank == 0){printf("Decomposition..........: OK \n");}

	buildLocalMesh(hLclData,nGlbElem,nGlbNodes,maxEl4GlbNd,h_nLclNodes,h_neig,h_neigFc,h_nFaces,h_part,h_boundCond,h_nEl4GlbNd,h_idGlbEl4GlbNd,h_elNd2GlbNd,h_cellCenter,h_fcCenter,h_glbNdCrdnts);
	if(rank == 0){printf("Build Local Mesh.......: OK \n");}

	//if(input.restartFlag==0){writeMeshRestart(hLclData, rank, size, nGlbElem, nGlbNodes);}
	//if(rank == 0){printf("Restart Mesh Files.....: OK \n");}

	startDevice(nGlbElem, nGlbNodes, rank, size, hLclData, mesh, field, comm, exp);
	if(rank == 0){printf("Device Started.........: OK \n");}

	if(input.restartFlag==0){
		initialize(mesh->nElem, mesh->cellCenter, field->consVec);
		iStart = 0;
		iStep = -1;
	}
	else if(input.restartFlag==1){
		readField(mesh->nElem, rank, field->consVec);
		iStep = iStart - 1;
	}
	if(rank == 0){printf("Initialize.............: OK \n");}

	writeFieldRestart(rank, mesh, field);
	if(rank == 0){printf("Restart Field Files....: OK \n");}

	export_plt(exp,field->consVec);
	if(rank == 0){printf("Export.................: OK \n");}
	
	size_t free_byte, total_byte;
    	cudaError_t status = cudaMemGetInfo(&free_byte, &total_byte);
    	if (status != cudaSuccess) {
        	printf("CUDA error: %s\n", cudaGetErrorString(status));
    	}

    	double free_gb = (double)free_byte / 1e9;
    	double total_gb = (double)total_byte / 1e9;
    	double used_gb = total_gb - free_gb;

    	if(rank==0){printf("GPU memory usage: %.2f GB used / %.2f GB total (%.2f GB free)\n",used_gb, total_gb, free_gb);}

}
