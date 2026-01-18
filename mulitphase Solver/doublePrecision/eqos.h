__device__ __forceinline__ double eqos(INPUT *d_input, int imat, int idVar3, double val1, double val2){

	double var3 = 0;
	switch(idVar3){ 
		case 0:
			var3 = (d_input->gamma[imat] - 1)*val1*val2 - d_input->gamma[imat]*d_input->pInf[imat] - (d_input->gamma[imat] - 1)*d_input->eta[imat]*val1;
		break;
		case 1:
			var3 = (d_input->gamma[imat] - 1)*val1*val2 - d_input->gamma[imat]*d_input->pInf[imat] - (d_input->gamma[imat] - 1)*d_input->eta[imat]*val1;
			var3 = (var3 + d_input->pInf[imat])/( (d_input->Cp[imat]/d_input->gamma[imat])*(d_input->gamma[imat] - 1)*val1 );
		break;
		case 2:
			var3 = (val1+d_input->pInf[imat])/((d_input->Cp[imat]/d_input->gamma[imat])*val2*(d_input->gamma[imat] - 1));
		break;
		case 3:
			var3 = (val1+d_input->gamma[imat]*d_input->pInf[imat] + (d_input->gamma[imat] - 1)*d_input->eta[imat]*val2)/(val2*(d_input->gamma[imat]-1));
		break;
		case 4:
			var3 = sqrt((d_input->gamma[imat]*(val2+d_input->pInf[imat]))/val1);
		break;
	}
	return var3;
}
