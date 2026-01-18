#include "strdata.h"

int main(int argc, char** argv){
    
	MESH   *mesh;
	FIELD  *field;
	EXPORT *exp;
	COMM   *comm;
	
	MPI_Init(&argc, &argv);
	
	mesh   = (MESH*)   malloc(sizeof(MESH));
	field  = (FIELD*)  malloc(sizeof(FIELD));
	exp    = (EXPORT*) malloc(sizeof(EXPORT));
	comm   = (COMM*)   malloc(sizeof(COMM));

	buildInput();
		
	mainStart(mesh, exp, comm, field);

	explicitLoop(mesh, exp, comm, field);	
    
	MPI_Finalize();

	exit(0);	
}

