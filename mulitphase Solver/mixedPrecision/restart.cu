#include "strdata.h"
/*
void writeMeshRestart(HOSTLOCALDATA *hLclData, int rank, int size, int nGlbElem, int nGlbNodes){
	
	struct stat st = {0};
	char pathRank[256];
	char pathMesh[256];
	char pathComm[256];
	char pathExp[256];
	char elementData[256];
	char faceData[256];
	char rankData[256];
	char nodeData[256];
	if(rank==0){
		if (stat("restartFiles", &st) == -1){
			mkdir("restartFiles", 0755);
		}
		if (stat("restartFiles/restartMesh", &st) == -1){
			mkdir("restartFiles/restartMesh", 0755);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(pathRank, "restartFiles/restartMesh/processor_%03d", rank);
	if (stat(pathRank, &st) == -1){
		mkdir(pathRank, 0755);
	}
	
//==================================== MESH STRUCT =========================================================//	
	sprintf(pathMesh, "%s/mesh", pathRank);
	if (stat(pathMesh, &st) == -1){
		mkdir(pathMesh, 0755);
	}

	sprintf(elementData, "%s/elementData.dat", pathMesh);
	FILE *f1 = fopen(elementData, "w");
	fprintf(f1, "%d %d %d\n",hLclData->nElem, hLclData->nGhostElem, hLclData->nOuter);
	for(int iel = 0; iel < hLclData->nElem; iel++){
		fprintf(f1, "%d %.10e %.10e %.10e %.10e\n", iel, hLclData->cellCenter[iel*3+0], hLclData->cellCenter[iel*3+1], hLclData->cellCenter[iel*3+2], hLclData->volume[iel]);
	}
	fclose(f1);

	sprintf(faceData, "%s/faceData.dat", pathMesh);
	FILE *f2 = fopen(faceData, "w");
	fprintf(f2, "%d %d\n",hLclData->nFcs, hLclData->nProxTot);
	for(int ifc = 0; ifc < hLclData->nFcs; ifc++){
		fprintf(f2, "%d %d %d %d %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", ifc,
				hLclData->fc2el[ifc*2+0],    hLclData->fc2el[ifc*2+1],    hLclData->boundCond[ifc], 
				hLclData->fcCenter[ifc*3+0], hLclData->fcCenter[ifc*3+1], hLclData->fcCenter[ifc*3+2], 
				hLclData->n[ifc*3+0],        hLclData->n[ifc*3+1],        hLclData->n[ifc*3+2], 
				hLclData->nt1[ifc*3+0],      hLclData->nt1[ifc*3+1],      hLclData->nt1[ifc*3+2],
				hLclData->nt2[ifc*3+0],      hLclData->nt2[ifc*3+1],      hLclData->nt2[ifc*3+2], 
				hLclData->area[ifc]);
	}
	fclose(f2);

//==================================== COMM STRUCT =========================================================//	
	sprintf(pathComm, "%s/comm", pathRank);
	if (stat(pathComm, &st) == -1){
		mkdir(pathComm, 0755);
	}
	sprintf(faceData, "%s/faceData.dat", pathComm);
	FILE *f3 = fopen(faceData, "w");
	fprintf(f3, "%d\n",hLclData->nProxFacesMax);
	for(int ifc = 0; ifc < hLclData->nFcs; ifc++){
		fprintf(f3,"%d %d %d %d\n", ifc, hLclData->neigRank4fc[ifc], hLclData->lclProx4fc[ifc], hLclData->lclFc2idRcv[ifc]); 
	}
	fclose(f3);

	sprintf(rankData, "%s/rankData.dat", pathComm);
	FILE *f4 = fopen(rankData, "w");
	fprintf(f4, "%d\n", hLclData->nNeigRanks);
	for(int irank = 0; irank < size; irank++){
		if(irank < hLclData->nNeigRanks){
			fprintf(f4, "%d %d %d\n", irank, hLclData->neigRanks[irank], hLclData->nProxFaces[irank]);
		}else{
			fprintf(f4, "%d %d %d\n", irank, -1, hLclData->nProxFaces[irank]);
		}
	}
	fclose(f4);

//==================================== EXP STRUCT =========================================================//	
	sprintf(pathExp, "%s/exp", pathRank);
	if (stat(pathExp, &st) == -1){
		mkdir(pathExp, 0755);
	}
	sprintf(elementData, "%s/elementData.dat", pathExp);
	FILE *f5 = fopen(elementData, "w");
	fprintf(f5, "%d %d\n",nGlbElem, hLclData->nElem);
	for(int iel = 0; iel < hLclData->nElem; iel++){
		fprintf(f5, "%d %d %d %d %d %d %d %d %d %d %d\n", iel, hLclData->lcl2glbEl[iel], 
			hLclData->elNd2lclNd[iel*8+0], hLclData->elNd2lclNd[iel*8+1], hLclData->elNd2lclNd[iel*8+2], hLclData->elNd2lclNd[iel*8+3], 
			hLclData->elNd2lclNd[iel*8+4], hLclData->elNd2lclNd[iel*8+5], hLclData->elNd2lclNd[iel*8+6], hLclData->elNd2lclNd[iel*8+7], 
			hLclData->nElNds[iel]); 
	}
	fclose(f5);

	sprintf(nodeData, "%s/nodeData.dat", pathExp);
	FILE *f6 = fopen(nodeData, "w");
	fprintf(f6, "%d %d %d\n", nGlbNodes, hLclData->nNodes, hLclData->maxEl4nd);
	for(int ind = 0; ind < hLclData->nNodes; ind++){
		fprintf(f6, "%d %d %d ", ind, hLclData->lcl2glbNd[ind], hLclData->nEl4nd[ind]);
		for(int id = 0; id < hLclData->maxEl4nd; id++){
			fprintf(f6, "%d ", hLclData->idEl4nd[ind*hLclData->maxEl4nd+id]);
		}
		fprintf(f6, "%.10e %.10e %.10e\n", hLclData->glbNdCord[ind*3+0], hLclData->glbNdCord[ind*3+1], hLclData->glbNdCord[ind*3+2]);
	}
	fclose(f6);

}
*/

void writeFieldRestart(int rank, MESH *mesh, FIELD *field){
	
	struct stat st = {0};
	char pathField[256];
	char pathRank[256];
	char file[256];
	
	double *h_consVec = (double*)malloc((mesh->nElem + mesh->nGhostElem)*input.NEQ*sizeof(double));
	cudaMemcpy(h_consVec, field->consVec, (mesh->nElem + mesh->nGhostElem)*input.NEQ*sizeof(double), cudaMemcpyDeviceToHost);

	sprintf(pathField, "restartFiles/restartField_%07d", iStep+1);
	if(rank==0){
		if (stat("restartFiles", &st) == -1){
			mkdir("restartFiles", 0755);
		}
		if (stat(pathField, &st) == -1){
			mkdir(pathField, 0755);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(pathRank, "%s/processor_%03d", pathField, rank);
	if (stat(pathRank, &st) == -1){
		mkdir(pathRank, 0755);
	}

	sprintf(file, "%s/consVecVal.dat", pathRank);
	FILE *f1 = fopen(file, "w");
	fprintf(f1, "%le\n", physicalTime);
	for(int iel = 0; iel < mesh->nElem; iel++){
		fprintf(f1, "%d",iel);
		for(int iVar = 0; iVar < input.NEQ; iVar++){
			fprintf(f1, " %.10e",h_consVec[iel*input.NEQ+iVar]);
		}
		fprintf(f1, "\n");
	}
	fclose(f1);

	free(h_consVec);	
}	
/*
void readMesh(HOSTLOCALDATA *hLclData, int rank){

	char pathRank[256];
	char pathMesh[256];
	char elementData[256];
	char faceData[256];

	sprintf(pathRank, "restartFiles/restartMesh/processor_%03d", rank);
	sprintf(pathMesh, "%s/mesh", pathRank);

	sprintf(elementData, "%s/elementData.dat", pathMesh);
	sprintf(faceData, "%s/faceData.dat", pathMesh);

	FILE *f1 = fopen(elementData, "r");
	if (!f1) {
		fprintf(stderr, "Error: cannot open %s\n", elementData);
		exit(EXIT_FAILURE);
	}

	// ================= ELEMENT DATA ==================
	fscanf(f1, "%d %d %d", &hLclData->nElem, &hLclData->nGhostElem, &hLclData->nOuter);

	// Allocate memory
	hLclData->cellCenter = (float *)malloc(hLclData->nElem * 3 * sizeof(float));
	hLclData->volume     = (float *)malloc((hLclData->nElem +hLclData->nGhostElem) *     sizeof(float));

	if (!hLclData->cellCenter || !hLclData->volume) {
		fprintf(stderr, "Memory allocation failed (elements)\n");
		exit(EXIT_FAILURE);
	}

	// Read element data lines
	for (int i = 0; i < hLclData->nElem; i++){
		int idx;
		float cx, cy, cz, vol;
		fscanf(f1, "%d %f %f %f %f", &idx, &cx, &cy, &cz, &vol);
			hLclData->cellCenter[idx * 3 + 0] = cx;
			hLclData->cellCenter[idx * 3 + 1] = cy;
			hLclData->cellCenter[idx * 3 + 2] = cz;
			hLclData->volume[idx] = vol;
	}
	fclose(f1);

	// ================= FACE DATA ==================
	FILE *f2 = fopen(faceData, "r");
	if (!f2) {
		fprintf(stderr, "Error: cannot open %s\n", faceData);
		exit(EXIT_FAILURE);
	}

	fscanf(f2, "%d %d", &hLclData->nFcs, &hLclData->nProxTot);

	// Allocate memory
	hLclData->fc2el     = (int *)malloc(hLclData->nFcs * 2 *   sizeof(int));
	hLclData->boundCond = (int *)malloc(hLclData->nFcs *       sizeof(int));
	hLclData->fcCenter  = (float *)malloc(hLclData->nFcs * 3 * sizeof(float));
	hLclData->n         = (float *)malloc(hLclData->nFcs * 3 * sizeof(float));
	hLclData->nt1       = (float *)malloc(hLclData->nFcs * 3 * sizeof(float));
	hLclData->nt2       = (float *)malloc(hLclData->nFcs * 3 * sizeof(float));
	hLclData->area      = (float *)malloc(hLclData->nFcs *     sizeof(float));

	if (!hLclData->fc2el || !hLclData->boundCond || !hLclData->fcCenter || 
	!hLclData->n || !hLclData->nt1 || !hLclData->nt2 || !hLclData->area) {
		fprintf(stderr, "Memory allocation failed (faces)\n");
		exit(EXIT_FAILURE);
	}

	// Read face data lines
	for (int i = 0; i < hLclData->nFcs; i++) {
		int ifc, el0, el1, bc;
		float fcx, fcy, fcz;
		float nx, ny, nz;
		float nt1x, nt1y, nt1z;
		float nt2x, nt2y, nt2z;
		float area;

		fscanf(f2, "%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f", 
			&ifc, &el0, &el1, &bc,
			&fcx, &fcy, &fcz,
			&nx, &ny, &nz,
			&nt1x, &nt1y, &nt1z,
			&nt2x, &nt2y, &nt2z,
			&area);

		hLclData->fc2el[ifc * 2 + 0] = el0;
		hLclData->fc2el[ifc * 2 + 1] = el1;
		hLclData->boundCond[ifc] = bc;

		hLclData->fcCenter[ifc * 3 + 0] = fcx;
		hLclData->fcCenter[ifc * 3 + 1] = fcy;
		hLclData->fcCenter[ifc * 3 + 2] = fcz;

		hLclData->n[ifc * 3 + 0] = nx;
		hLclData->n[ifc * 3 + 1] = ny;
		hLclData->n[ifc * 3 + 2] = nz;

		hLclData->nt1[ifc * 3 + 0] = nt1x;
		hLclData->nt1[ifc * 3 + 1] = nt1y;
		hLclData->nt1[ifc * 3 + 2] = nt1z;
	
		hLclData->nt2[ifc * 3 + 0] = nt2x;
		hLclData->nt2[ifc * 3 + 1] = nt2y;
		hLclData->nt2[ifc * 3 + 2] = nt2z;

		hLclData->area[ifc] = area;
	}

	fclose(f2);
}

void readComm(HOSTLOCALDATA *hLclData, int rank, int size){

	char pathRank[256];
	char pathComm[256];
	char faceData[256];
	char rankData[256];

	sprintf(pathRank, "restartFiles/restartMesh/processor_%03d", rank);
	sprintf(pathComm, "%s/comm", pathRank);

	sprintf(faceData, "%s/faceData.dat", pathComm);
	sprintf(rankData, "%s/rankData.dat", pathComm);

	FILE *f3 = fopen(faceData, "r");
	if (!f3) {
		fprintf(stderr, "Error: cannot open %s\n", faceData);
		exit(EXIT_FAILURE);
	}

	// ---------- FACE DATA ----------
	fscanf(f3, "%d", &hLclData->nProxFacesMax);

	// Allocate arrays
	hLclData->neigRank4fc = (int *)malloc(hLclData->nFcs * sizeof(int));
	hLclData->lclProx4fc  = (int *)malloc(hLclData->nFcs * sizeof(int));
	hLclData->lclFc2idRcv = (int *)malloc(hLclData->nFcs * sizeof(int));

	if (!hLclData->neigRank4fc || !hLclData->lclProx4fc || !hLclData->lclFc2idRcv){
		fprintf(stderr, "Memory allocation failed (comm face arrays)\n");
		exit(EXIT_FAILURE);
	}

	// Read each face line
	for (int i = 0; i < hLclData->nFcs; i++) {
		int ifc, neigR, prox, idRcv;
		fscanf(f3, "%d %d %d %d", &ifc, &neigR, &prox, &idRcv);
			hLclData->neigRank4fc[ifc] = neigR;
			hLclData->lclProx4fc[ifc]  = prox;
			hLclData->lclFc2idRcv[ifc] = idRcv;
	}
	
	fclose(f3);

	// ---------- RANK DATA ----------
	FILE *f4 = fopen(rankData, "r");
	if (!f4) {
		fprintf(stderr, "Error: cannot open %s\n", rankData);
		exit(EXIT_FAILURE);
	}

	fscanf(f4, "%d", &hLclData->nNeigRanks);

	// Allocate arrays
	hLclData->neigRanks  = (int *)malloc(hLclData->nNeigRanks * sizeof(int));
	hLclData->nProxFaces = (int *)malloc(size * sizeof(int));

	if (!hLclData->neigRanks || !hLclData->nProxFaces) {
		fprintf(stderr, "Memory allocation failed (comm rank arrays)\n");
		exit(EXIT_FAILURE);
	}

	for (int irank = 0; irank < size; irank++) {
		int id, neig, nprox;
		fscanf(f4, "%d %d %d", &id, &neig, &nprox);
		hLclData->nProxFaces[irank] = nprox;
		if (irank < hLclData->nNeigRanks){
			hLclData->neigRanks[irank] = neig;
		}
	}

	fclose(f4);
}

void readExp(HOSTLOCALDATA *hLclData, int rank, int *nGlbNodes, int *nGlbElem){

	char pathRank[256];
	char pathExp[256];
	char elementData[256];
	char nodeData[256];

	sprintf(pathRank, "restartFiles/restartMesh/processor_%03d", rank);
	sprintf(pathExp, "%s/exp", pathRank);

	sprintf(elementData, "%s/elementData.dat", pathExp);
	sprintf(nodeData, "%s/nodeData.dat", pathExp);

	// ====================== ELEMENT DATA ==========================
	FILE *f5 = fopen(elementData, "r");
	if (!f5) {
		fprintf(stderr, "Error: cannot open %s\n", elementData);
		exit(EXIT_FAILURE);
	}

	fscanf(f5, "%d %d", nGlbElem, &hLclData->nElem);

	hLclData->lcl2glbEl  = (int *)malloc(hLclData->nElem * sizeof(int));
	hLclData->nElNds     = (int *)malloc(hLclData->nElem * sizeof(int));
	hLclData->elNd2lclNd = (int *)malloc(hLclData->nElem * 8 * sizeof(int));

	if (!hLclData->lcl2glbEl || !hLclData->nElNds || !hLclData->elNd2lclNd) {
		fprintf(stderr, "Memory allocation failed (EXP element arrays)\n");
		exit(EXIT_FAILURE);
	}

	for (int iel = 0; iel < hLclData->nElem; iel++) {
		int idx, lcl2glb, nd[8], nnds;
		fscanf(f5, "%d %d %d %d %d %d %d %d %d %d %d",
			&idx, &lcl2glb,
			&nd[0], &nd[1], &nd[2], &nd[3],
			&nd[4], &nd[5], &nd[6], &nd[7],
			&nnds);

		hLclData->lcl2glbEl[idx] = lcl2glb;
		for (int j = 0; j < 8; j++){
			hLclData->elNd2lclNd[idx * 8 + j] = nd[j];
		}
		hLclData->nElNds[idx] = nnds;
	}
	fclose(f5);

	// ====================== NODE DATA ==========================
	FILE *f6 = fopen(nodeData, "r");
	if (!f6) {
		fprintf(stderr, "Error: cannot open %s\n", nodeData);
		exit(EXIT_FAILURE);
	}

	fscanf(f6, "%d %d %d", nGlbNodes, &hLclData->nNodes, &hLclData->maxEl4nd);

	hLclData->lcl2glbNd = (int *)malloc(hLclData->nNodes * sizeof(int));
	hLclData->nEl4nd    = (int *)malloc(hLclData->nNodes * sizeof(int));
	hLclData->idEl4nd   = (int *)malloc(hLclData->nNodes * hLclData->maxEl4nd * sizeof(int));
	hLclData->glbNdCord = (float *)malloc(hLclData->nNodes * 3 * sizeof(float));

	if (!hLclData->lcl2glbNd || !hLclData->nEl4nd || !hLclData->idEl4nd || !hLclData->glbNdCord) {
		fprintf(stderr, "Memory allocation failed (EXP node arrays)\n");
		exit(EXIT_FAILURE);
	}

	for (int ind = 0; ind < hLclData->nNodes; ind++) {
		int idx, l2g, nels;
		fscanf(f6, "%d %d %d", &idx, &l2g, &nels);
		hLclData->lcl2glbNd[idx] = l2g;
		hLclData->nEl4nd[idx] = nels;

		for (int id = 0; id < hLclData->maxEl4nd; id++) {
			int val;
			fscanf(f6, "%d", &val);
			hLclData->idEl4nd[idx * hLclData->maxEl4nd + id] = val;
		}

		float x, y, z;
		fscanf(f6, "%f %f %f", &x, &y, &z);
		hLclData->glbNdCord[idx * 3 + 0] = x;
		hLclData->glbNdCord[idx * 3 + 1] = y;
		hLclData->glbNdCord[idx * 3 + 2] = z;
	}	
	fclose(f6);
}
*/
void readField(int nElem, int rank, double *d_consVec){
	
	double *h_consVec = (double*)malloc(nElem*input.NEQ*sizeof(double));

	// 1. find the latest restartField directory by scanning "restartFiles"
	DIR *d = opendir("restartFiles");
	if (!d) {
		fprintf(stderr, "Error: cannot open restartFiles directory\n");
		exit(EXIT_FAILURE);
	}
	struct dirent *entry;
	int bestStep = -1;
	char bestDirName[256] = {0};

	while ((entry = readdir(d)) != NULL){
		const char *prefix = "restartField_";
		if (strncmp(entry->d_name, prefix, strlen(prefix)) == 0){
			int step = atoi(entry->d_name + strlen(prefix));
			if (step > bestStep) {
				bestStep = step;
				strncpy(bestDirName, entry->d_name, sizeof(bestDirName)-1);
			}
		}
	}
	closedir(d);

	if (bestStep < 0) {
		fprintf(stderr, "No restartField_xxxxx directories found\n");
		exit(EXIT_FAILURE);
	}

	if(rank == 0){printf("Restart Step...........: %d\n", bestStep);}
	iStart = bestStep;
	
	// 2. build path to the consVec file
	char pathRank[256];
	char filepath[512];
	snprintf(pathRank, sizeof(pathRank), "restartFiles/%s/processor_%03d",bestDirName, rank);
	snprintf(filepath, sizeof(filepath), "%s/consVecVal.dat", pathRank);

	// 3. open file and read lines
	FILE *f = fopen(filepath, "r");
	if (!f) {
		fprintf(stderr, "Error: cannot open %s\n", filepath);
		exit(EXIT_FAILURE);
	}

	// 4. read data
	fscanf(f, "%lf",&physicalTime);
	for (int iel = 0; iel < nElem; iel++) {
		int idx;
		double v0, v1, v2, v3, v4, v5, v6, v7, v8, v9;
		int ret = fscanf(f, "%d %lf %lf %lf %lf %lf",&idx, &v0, &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9);
		if (ret != 6) {
			fprintf(stderr, "Error reading line %d in %s\n", iel, filepath);
			exit(EXIT_FAILURE);
		}
		h_consVec[idx * 5 + 0] = v0;
		h_consVec[idx * 5 + 1] = v1;
		h_consVec[idx * 5 + 2] = v2;
		h_consVec[idx * 5 + 3] = v3;
		h_consVec[idx * 5 + 4] = v4;
		h_consVec[idx * 5 + 5] = v5;
		h_consVec[idx * 5 + 6] = v6;
		h_consVec[idx * 5 + 7] = v7;
		h_consVec[idx * 5 + 8] = v8;
		h_consVec[idx * 5 + 9] = v9;
	}
	fclose(f);
	
	cudaMemcpy(d_consVec, h_consVec,nElem*10*sizeof(double), cudaMemcpyHostToDevice);
	free(h_consVec);

}

void restartDomain(MESH *mesh, EXPORT *exp, COMM *comm, FIELD *field){
/*	
	int nGlbNodes, nGlbElem;
	
	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size); 

	HOSTLOCALDATA *hLclData   = (HOSTLOCALDATA*) malloc(sizeof(HOSTLOCALDATA));

	readMesh(hLclData, rank);
	if(rank == 0){printf("Restart Mesh Struct....: OK \n");}

	readComm(hLclData, rank, size);
	if(rank == 0){printf("Restart Comm Struct....: OK \n");}

	readExp(hLclData, rank, &nGlbNodes, &nGlbElem);
	if(rank == 0){printf("Restart Exp  Struct....: OK \n");}
	
	communicateVolumes(hLclData, rank, size, hLclData->nNeigRanks, hLclData->nProxFacesMax, hLclData->neigRanks, hLclData->nProxFaces, hLclData->neigRank4fc, hLclData->lclProx4fc, hLclData->lclFc2idRcv, hLclData->fc2el, hLclData->volume);
			
	startDevice(nGlbElem, nGlbNodes, rank, size, hLclData, mesh, field, comm, exp);
	if(rank == 0){printf("Device Started.........: OK \n");}
	
	readField(mesh->nElem, rank, field->consVec);
	if(rank == 0){printf("Restart Field..........: OK \n");}
		
	size_t free_byte, total_byte;
    	cudaError_t status = cudaMemGetInfo(&free_byte, &total_byte);
    	if (status != cudaSuccess) {
        	printf("CUDA error: %s\n", cudaGetErrorString(status));
    	}

    	float free_gb = (float)free_byte / 1e9;
    	float total_gb = (float)total_byte / 1e9;
    	float used_gb = total_gb - free_gb;

    	if(rank==0){printf("GPU memory usage: %.2f GB used / %.2f GB total (%.2f GB free)\n",used_gb, total_gb, free_gb);}	
	MPI_Barrier(MPI_COMM_WORLD);
*/	
}
