#include "strdata.h"

INPUT *d_input;  
INPUT input;     
int iStep = 0;
int iStart = 0;
double physicalTime = 0.0;
int printHLLC=0;

void buildInput(){
	
	input.dt	 = 0.0;//3.7e-9;
	input.adaptDt    = 1;
	input.adaptDtItr = 100; 
	input.maxDt	 = 1e-8;//1.00e-9; worked
	input.cfl	 = 0.1;//0.1;
	
	input.avfMin	 = 1e-4;
	
	input.rkSteps		= 1;
	input.finalTime		= 3.5e-3;
	input.exportTime	= 2.0e-5;
	input.restartStep	= 100e3;
	input.restartFlag	= 0;
	input.neglectX		= 0;
	input.neglectY		= 0;
	input.neglectZ		= 0;
	input.nMaterials	= 2;
	input.NEQ		= 10;
	input.scale		= 0.001;
	input.viscFlag		= 1; 
	input.turbulence	= 1;
	input.woodsFlag		= 1;
	
	input.thermFntRlxFlag	= 0;
	input.thermRlxRate	= 0.5;
	input.thermUrf		= 0.8;	
	input.thermItrMax	= 25;   	

	input.chemFlag		= 1;
	input.chemFntRlxFlag	= 0;
	input.chemRlxRate	= 0;
	input.chemUrf		= 0.5;
	input.chemItrMax	= 40;
	input.chemCondesation	= 0;
	input.chemLimit		= 20e-4;

	input.consCalcFlag	= 1;
	input.consCalcIvar	= 9;

	input.satTempUrf	= 0.8;
	
	input.gamma[0]		= 2.17290;
	input.pInf[0]		= 114907765; 
	input.Cp[0]		= 1854.4;       
	input.eta[0]		= -302568.5557; 
	input.satConst[0]	= -488.9933492;
	input.dynVisc[0]	= 1.8e-4;      
	
	input.gamma[1]		= 1.4588;
	input.pInf[1]		= 0.0;
	input.Cp[1]		= 965.98;
	input.eta[1]		= 9741.04362;	
	input.satConst[1]	= 297.63826;
	input.dynVisc[1]	= 1.99e-5; 	

	cudaMalloc((void**)&d_input, sizeof(INPUT));
	cudaMemcpy(d_input, &input, sizeof(INPUT), cudaMemcpyHostToDevice);
}

