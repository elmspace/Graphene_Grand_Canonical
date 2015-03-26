void omega(double ****w){

  int i,j,k;
  int ii,jj,kk;
  double junk;
 

  if(Iomega==0){
    
    std::ifstream infile;
    if(AlphaB==1){infile.open("./OMEGA/omega_20_20_20_phiC_99.read");}
    if(Bilayer==1){infile.open("./OMEGA/omega_20_20_20_Bilayer.read");}
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  infile >> ii >> jj >> kk >> w[0][i][j][k] >> w[1][i][j][k] >> w[2][i][j][k] >> w[3][i][j][k] >> w[4][i][j][k] >> w[5][i][j][k];
	}
      }
    }
    infile.close();

  }else if(Iomega==1){
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ Make Structure
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){

	  w[0][i][j][k]=-10.0*cos(i*2.0*Pi/Nx); //A
	  w[1][i][j][k]=10.0*cos(i*2.0*Pi/Nx); //C
	  w[2][i][j][k]=0.0; //B1
	  w[3][i][j][k]=0.0; //B2
	  w[4][i][j][k]=0.0; //B3
	  w[5][i][j][k]=0.0; //B4
	  
	}
      }
    }
 
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  }else if(Iomega==2){
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){

	  w[0][i][j][k]=-50.0*(drand48()-0.50); //A1
	  w[1][i][j][k]=-50.0*(drand48()-0.50); //A2
	  w[2][i][j][k]=-50.0*(drand48()-0.50); //B1
	  w[3][i][j][k]=-50.0*(drand48()-0.50); //B2

	}
      }
    }
  


  }
 
};



