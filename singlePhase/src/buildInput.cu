#include "strdata.h"

INPUT *d_input;  
INPUT input;     
int iStep = 0;
int iStart = 0;
float physicalTime = 0.0;
int printHLLC=0;

void buildInput(){
	
	input.rkSteps     = 1; 
	input.nSteps      = 2e6;
	input.exportStep  = 25e3;
	input.restartStep = 250e3;
	input.restartFlag = 0;
	input.neglectX    = 0;
	input.neglectY    = 0;
	input.neglectZ    = 0;
	input.nMaterials  = 1;
	input.NEQ         = 5;
	input.scale       = 0.001;
	input.turbulence  = 1;

	input.dt          = 5e-9;
	input.gamma       = 1.4;
	input.pInf        = 0.0;
	input.Cp          = 1006;
	input.dynVisc     = 1.9212e-5;

	cudaMalloc((void**)&d_input, sizeof(INPUT));
	cudaMemcpy(d_input, &input, sizeof(INPUT), cudaMemcpyHostToDevice);
}

