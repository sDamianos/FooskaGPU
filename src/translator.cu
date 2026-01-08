#include "strdata.h"

void translator_ptw(float scale, int *nGlbElem, int *nGlbNodes, int **nFaces, int **nLclNodes, int **elNd2GlbNd, int **boundCond, float **cellCenter, float **fcCenter, float **glbNdCrdnts) {

	int nd,el,typ,numcon,c1,c2,c3,c4,c5,c6,c7,c8;
	int und1,und2,und3,bounds,topic_face,bccon;
	int ndnum,elnum,bcnum,dim,code1,code2;
	float xc,yc,zc;
	int i,j,x;
	char line[100], text[50];
	FILE *neu;


	neu    = fopen("sh.neu","r");

	if (neu==NULL) {
		printf("ForestFV: TRANSLATOR: mesh file not found exiting \n ");
		exit(0);
	}
  
	memset(line, 0, sizeof(line));
 

	char CN_INFO[100]={'C','O','N','T','R','O','L',' ','I','N','F','O',' ',' ',' ',' ',' ',' ',' ',' ',' ','1','.','2','.','1' ,'\r','\n'};
	char ND_COOR[100]={'N','O','D','A','L',' ','C','O','O','R','D','I','N','A','T','E','S',' ',' ',' ',' ','1','.','2','.','1','\r','\n'};
	char EL_CONN[100]={'E','L','E','M','E','N','T','S','/','C','E','L','L','S',' ',' ',' ',' ',' ',' ',' ','1','.','2','.','1','\r','\n'};
	char BC_TEX[100]={'B','O','U','N','D','A','R','Y',' ','C','O','N','D','I','T','I','O','N','S',' ','1','.','2','.','1','\r','\n'};
  

//=====================================================================================================
	//  Header
	j=0;
	while((x=strncmp(CN_INFO,line,20))!=0){   // Loops until "CONTROL INFO         1.2.1"
		for(i=0;i<100;i++){
			line[i]=0;
		}
		fgets(line,100,neu);
		if (j>10000) {
			printf("1: Wrong translator for current mesh. Gambit mesh: 1 Pointwise: 2 \n");
			exit(0);
		}
		j++; 
	}; 

	fgets(line,100,neu);
	fgets(line,100,neu);
	fgets(line,100,neu);
	fgets(line,100,neu);
	fgets(line,100,neu);

	fscanf(neu,"%d \t %d \t %d \t %d \t %d \t %d",&ndnum,&elnum,&code1,&bcnum,&dim,&code2);


	j=0;
	while((x=strncmp(ND_COOR,line,20))!=0){   // Loops until "NODAL COORDINATES    1.2.1"
		for(i=0;i<100;i++){
			line[i]=0;
		}
		fgets(line,100,neu); 
		if (j>10000) {
			printf("2: Wrong translator for current mesh. Gambit mesh: 1 Pointwise: 2 \n");
			exit(0);
		}
		j++; 
	};

//=====================================================================================================
// Nodes
  
	//Construction of .nd file        
	(*glbNdCrdnts)  = (float*) malloc((ndnum*3)*sizeof(float));
	
	for(i=0;i<ndnum;i++){
		fscanf(neu,"%d\t%f\t%f\t%f\n",&nd,&xc,&yc,&zc);
		(*glbNdCrdnts)[i*3+0] = xc*scale; 
		(*glbNdCrdnts)[i*3+1] = yc*scale;
		(*glbNdCrdnts)[i*3+2] = zc*scale;
	}
  
	j=0;
	while((x=strncmp(EL_CONN,line,20))!=0){  // Loops until "NODAL COORDINATES    1.2.1"
		for(i=0;i<100;i++){
			line[i]=0;
		}
		fgets(line,100,neu); 
		if (j>10000) {
			printf("3: Wrong translator for current mesh. Gambit mesh: 1 Pointwise: 2 \n");
			exit(0);
		}
		j++; 
	};

//=====================================================================================================
// Elements/Cells, 20012
  
  	
	*cellCenter = (float*) malloc((elnum*3)*sizeof(float));
	*fcCenter   = (float*) malloc((elnum*6*3)*sizeof(float));
	*elNd2GlbNd  = (int*) malloc((elnum*8)*sizeof(int));
	*boundCond  = (int*) malloc((elnum*6)*sizeof(int));
	*nFaces     = (int*) malloc(elnum*sizeof(int));
	*nLclNodes  = (int*) malloc(elnum*sizeof(int));
	*nGlbNodes = ndnum;
	*nGlbElem = elnum;
	for(i=0; i<elnum; i++){
		(*nFaces)[i]     = 0;
		(*nLclNodes)[i]  = 0;
		(*cellCenter)[i] = 0.0;
		for(j=0; j<8; j++){
			(*elNd2GlbNd)[i*8+j] = 0;
			if(j<6){
				(*boundCond)[i*6+j] = 0;
				for(int iv = 0; iv < 3; iv++){
					(*fcCenter)[i*6+j*3+iv]  = 0.0;
				}
			}if(j<3){
				(*cellCenter)[i*3+j] = 0.0;
			}
		}	
	}


	for(i=0;i<elnum;i++){
		fscanf(neu,"%d\t%d\t%d",&el,&typ,&numcon);
		
		if(numcon==8){  //Linear Hex
			fscanf(neu,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8);
			(*elNd2GlbNd)[(el-1)*8 + 0] = c1-1;
			(*elNd2GlbNd)[(el-1)*8 + 1] = c2-1;
			(*elNd2GlbNd)[(el-1)*8 + 2] = c3-1;
			(*elNd2GlbNd)[(el-1)*8 + 3] = c4-1;
			(*elNd2GlbNd)[(el-1)*8 + 4] = c5-1;
			(*elNd2GlbNd)[(el-1)*8 + 5] = c6-1;
			(*elNd2GlbNd)[(el-1)*8 + 6] = c7-1;
			(*elNd2GlbNd)[(el-1)*8 + 7] = c8-1;
			(*nFaces)[el-1] = 6;
			(*nLclNodes)[el-1] = 8;
			
			(*cellCenter)[ (el-1)*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] + 
							  (*glbNdCrdnts)[(c5-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] + (*glbNdCrdnts)[(c7-1)*3 + 0] + (*glbNdCrdnts)[(c8-1)*3 + 0] )/8;
			(*cellCenter)[ (el-1)*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] + 
							  (*glbNdCrdnts)[(c5-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] + (*glbNdCrdnts)[(c7-1)*3 + 1] + (*glbNdCrdnts)[(c8-1)*3 + 1] )/8;
			(*cellCenter)[ (el-1)*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] + 
							  (*glbNdCrdnts)[(c5-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] + (*glbNdCrdnts)[(c7-1)*3 + 2] + (*glbNdCrdnts)[(c8-1)*3 + 2] )/8;

				     //Cell     ifc  xyz 
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 1*3 + 0 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c8-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 1 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c8-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 2 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c8-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] )/4;
			
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 0 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c7-1)*3 + 0] + (*glbNdCrdnts)[(c8-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 1 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c7-1)*3 + 1] + (*glbNdCrdnts)[(c8-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 2 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c7-1)*3 + 2] + (*glbNdCrdnts)[(c8-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 3*3 + 0 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] + (*glbNdCrdnts)[(c7-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 1 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] + (*glbNdCrdnts)[(c7-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 2 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] + (*glbNdCrdnts)[(c7-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 4*3 + 0 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 4*3 + 1 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 4*3 + 2 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 5*3 + 0 ] = ( (*glbNdCrdnts)[(c5-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] + (*glbNdCrdnts)[(c8-1)*3 + 0] + (*glbNdCrdnts)[(c7-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 5*3 + 1 ] = ( (*glbNdCrdnts)[(c5-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] + (*glbNdCrdnts)[(c8-1)*3 + 1] + (*glbNdCrdnts)[(c7-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 5*3 + 2 ] = ( (*glbNdCrdnts)[(c5-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] + (*glbNdCrdnts)[(c8-1)*3 + 2] + (*glbNdCrdnts)[(c7-1)*3 + 2] )/4;

		}

		if(numcon==6){  //Linear Prism
			fscanf(neu,"%d\t%d\t%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4,&c5,&c6);
			(*elNd2GlbNd)[(el-1)*8 + 0] = c1-1;
			(*elNd2GlbNd)[(el-1)*8 + 1] = c2-1;
			(*elNd2GlbNd)[(el-1)*8 + 2] = c3-1;
			(*elNd2GlbNd)[(el-1)*8 + 3] = c4-1;
			(*elNd2GlbNd)[(el-1)*8 + 4] = c5-1;
			(*elNd2GlbNd)[(el-1)*8 + 5] = c6-1;
			(*nFaces)[el-1] = 5;
			(*nLclNodes)[el-1] = 6;

			(*cellCenter)[ (el-1)*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + 
							  (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] )/6;
			(*cellCenter)[ (el-1)*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + 
							  (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] )/6;
			(*cellCenter)[ (el-1)*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + 
							  (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] )/6;
				     //Cell     ifc  xyz 
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 1*3 + 0 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 1 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 2 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 2*3 + 0 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] )/4;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 1 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] )/4;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 2 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 3*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] )/3; 
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] )/3; 
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] )/3; 

			(*fcCenter)[ (el-1)*6*3 + 4*3 + 0 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] + (*glbNdCrdnts)[(c6-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 4*3 + 1 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] + (*glbNdCrdnts)[(c6-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 4*3 + 2 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] + (*glbNdCrdnts)[(c6-1)*3 + 2] )/3;
			
		}

		if(numcon==5){  //Linear Pyramid
			fscanf(neu,"%d\t%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4,&c5);
			(*elNd2GlbNd)[(el-1)*8 + 0] = c1-1;
			(*elNd2GlbNd)[(el-1)*8 + 1] = c2-1;
			(*elNd2GlbNd)[(el-1)*8 + 2] = c3-1;
			(*elNd2GlbNd)[(el-1)*8 + 3] = c4-1;
			(*elNd2GlbNd)[(el-1)*8 + 4] = c5-1;
			(*nFaces)[el-1] = 5;
			(*nLclNodes)[el-1] = 5;

			(*cellCenter)[ (el-1)*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] +
							  (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] )/5;
			(*cellCenter)[ (el-1)*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] +
							  (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] )/5;
			(*cellCenter)[ (el-1)*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] +
							  (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] )/5;
				     //Cell     ifc  xyz 
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0]) /4;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1]) /4;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2]) /4;

			(*fcCenter)[ (el-1)*6*3 + 1*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0])/3;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1])/3;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2])/3;

			(*fcCenter)[ (el-1)*6*3 + 2*3 + 0 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 1 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 2 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] )/3;
			
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 0 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 1 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 2 ] = ( (*glbNdCrdnts)[(c4-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] )/3;

			(*fcCenter)[ (el-1)*6*3 + 4*3 + 0 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c5-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 4*3 + 1 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c5-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 4*3 + 2 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c5-1)*3 + 2] )/3;

		}	

		if(numcon==4){  //Linear Tetrahedral
			fscanf(neu,"%d\t%d\t%d\t%d\n",&c1,&c2,&c3,&c4);
			(*elNd2GlbNd)[(el-1)*8 + 0] = c1-1;
			(*elNd2GlbNd)[(el-1)*8 + 1] = c2-1;
			(*elNd2GlbNd)[(el-1)*8 + 2] = c3-1;
			(*elNd2GlbNd)[(el-1)*8 + 3] = c4-1;
			(*nFaces)[el-1] = 4;
			(*nLclNodes)[el-1] = 4;

			(*cellCenter)[ (el-1)*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] )/4;
			(*cellCenter)[ (el-1)*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] )/4;
			(*cellCenter)[ (el-1)*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] )/4;

			(*fcCenter)[ (el-1)*6*3 + 0*3 + 0 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 1 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 0*3 + 2 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] )/3;
			
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 0 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 1 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 1*3 + 2 ] = ( (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] )/3;
			
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 0 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 0] + (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 1 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 1] + (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 2*3 + 2 ] = ( (*glbNdCrdnts)[(c2-1)*3 + 2] + (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] )/3;

			(*fcCenter)[ (el-1)*6*3 + 3*3 + 0 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 0] + (*glbNdCrdnts)[(c1-1)*3 + 0] + (*glbNdCrdnts)[(c4-1)*3 + 0] )/3;
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 1 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 1] + (*glbNdCrdnts)[(c1-1)*3 + 1] + (*glbNdCrdnts)[(c4-1)*3 + 1] )/3;
			(*fcCenter)[ (el-1)*6*3 + 3*3 + 2 ] = ( (*glbNdCrdnts)[(c3-1)*3 + 2] + (*glbNdCrdnts)[(c1-1)*3 + 2] + (*glbNdCrdnts)[(c4-1)*3 + 2] )/3;
			

		}

	}

//=====================================================================================================

	// Skip to BC
	while((x=strncmp(BC_TEX,line,20))!=0){
		for(i=0;i<100;i++){
			line[i]=0;
		}
		fgets(line,100,neu);
	};


//=====================================================================================================
	//Boundary conditions
	int bc_symm = 0; int bc_out = 0; int bc_inl =0; int bc_slip = 0; int bc_wall = 0;
	for(bccon=0;bccon<bcnum;bccon++){

		fscanf(neu,"%s\t%d\t%d\t%d\t%d\n",text,&und1,&bounds,&und2,&und3);

		if     (und3==45){ bc_symm  = 1;} //symmetry file
		else if(und3==26){ bc_out   = 1;} // outlet
		else if(und3==13){ bc_inl   = 1;} //inlet
		else if(und3==42){ bc_slip  = 1;} //slip wall
		else if(und3==51){ bc_wall  = 1;} //wall

		if(und3==45){
			for(i=0;i<bounds;i++){                  
				fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
				(*boundCond)[(el-1)*6+(topic_face-1)] = 1;
			}
		}

		if(und3==26){
			for(i=0;i<bounds;i++){
				fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
				(*boundCond)[(el-1)*6+(topic_face-1)] = 2;
			}
		}

		if(und3==13){
			for(i=0;i<bounds;i++){
				fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
				(*boundCond)[(el-1)*6+(topic_face-1)] = 3;
			}
		}

		if(und3==42){
			for(i=0;i<bounds;i++){
				fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
				(*boundCond)[(el-1)*6+(topic_face-1)] = 4;
			}
		}

		if(und3==51){
			for(i=0;i<bounds;i++){
				fscanf(neu,"%d\t%d\t%d\n",&el,&typ,&topic_face);
				(*boundCond)[(el-1)*6+(topic_face-1)] = 5;
			}
		}

		if(bccon!=bcnum-1){
			fgets(line,100,neu);
			fgets(line,100,neu);
		}
	}//bcon

	fclose(neu);  
	return;

}
