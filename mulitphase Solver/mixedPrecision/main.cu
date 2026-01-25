#include "strdata.h"
void assign_device_to_rank() {

	int global_rank, local_rank, num_devices;
	MPI_Comm node_comm;

	MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, global_rank, MPI_INFO_NULL, &node_comm);

	MPI_Comm_rank(node_comm, &local_rank);

	cudaGetDeviceCount(&num_devices);
	cudaSetDevice(local_rank % num_devices);

	MPI_Comm_free(&node_comm);
}

int main(int argc, char** argv){
    
	MESH   *mesh;
	FIELD  *field;
	EXPORT *exp;
	COMM   *comm;
	
	MPI_Init(&argc, &argv);

	assign_device_to_rank();

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

