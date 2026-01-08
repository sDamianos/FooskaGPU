
__device__ __forceinline__ void rotateVector_glb2lcl(float vecGlb[3], float vecLcl[3], float * n, float * nt1, float * nt2){

	vecLcl[0] = vecGlb[0]*n[0] + vecGlb[1]*n[1] + vecGlb[2]*n[2];
	vecLcl[1] = vecGlb[0]*nt1[0] + vecGlb[1]*nt1[1] + vecGlb[2]*nt1[2]; 
	vecLcl[2] = vecGlb[0]*nt2[0] + vecGlb[1]*nt2[1] + vecGlb[2]*nt2[2];
}

__device__ __forceinline__ void rotateVector_lcl2glb(float vecLcl[3], float vecGlb[3], float * n, float * nt1, float * nt2){

	vecGlb[0] = vecLcl[0]*n[0] + vecLcl[1]*nt1[0] + vecLcl[2]*nt2[0];
	vecGlb[1] = vecLcl[0]*n[1] + vecLcl[1]*nt1[1] + vecLcl[2]*nt2[1]; 
	vecGlb[2] = vecLcl[0]*n[2] + vecLcl[1]*nt1[2] + vecLcl[2]*nt2[2];
}
