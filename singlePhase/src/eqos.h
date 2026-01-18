
__device__ __forceinline__ float eqos(INPUT *d_input, int idVar3, float val1, float val2){

	float var3 = 0;
	switch(idVar3){ 
		case 0:
			var3 = (d_input->gamma - 1)*val1*val2 - d_input->gamma*d_input->pInf;
		break;
		case 1:
			var3 = (d_input->gamma/d_input->Cp)*(val2 - d_input->pInf/val1);
		break;
		case 2:
			var3 = (val1+d_input->pInf)/((d_input->Cp/d_input->gamma)*val2*(d_input->gamma - 1));
		break;
		case 3:
			var3 = (val1+d_input->gamma*d_input->pInf)/(val2*(d_input->gamma-1));
		break;
		case 4:
			var3 = sqrt((d_input->gamma*(val2+d_input->pInf))/val1);
		break;
	}
	return var3;
}
