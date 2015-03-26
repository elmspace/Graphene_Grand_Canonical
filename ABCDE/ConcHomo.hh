double ConcHomo(double ****phi,double ****w,double *Ns,double ds,double ***k_vector,double *dxyz){

  int         i,j,l,s;
  double      Q;

  double      ****qB4;
  double      ***qintB4;


  qB4=create_4d_double_array(Nx,Ny,Nz,((int)Ns[5]+1),"qB4");

  qintB4=create_3d_double_array(Nx,Ny,Nz,"qintB4");

  //+++++++++++++++++++++++++++++++++++++++Forward++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB4[i][j][l]=1.0;
      }
    }
  }
  solveModDiffEqn_FFT(qB4,w[5],qintB4,ds,(int)Ns[5],1,k_vector,dxyz);
 

  //++++++++++++++++++++++++++++++++++++++Single Chain Partition Function+++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Q=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	Q+=(qB4[i][j][l][(int)Ns[5]])*dxyz[0]*dxyz[1]*dxyz[2];
      }
    }
  }
  // Normalizing with respect to the volume of the box
  Q/=((dxyz[0]*Nx)*(dxyz[1]*Ny)*(dxyz[2]*Nz));

  
  // Here we do the concentration calculation
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){

	phi[5][i][j][l]=0.0;  //B4

	//B4
	for(s=0;s<(Ns[5]+1);s++){
	  if(s==0 || s==(int)Ns[5]){
	    phi[5][i][j][l]+=0.5*qB4[i][j][l][s]*qB4[i][j][l][(int)Ns[5]-s]*ds;
	  }else{
	    phi[5][i][j][l]+=qB4[i][j][l][s]*qB4[i][j][l][(int)Ns[5]-s]*ds;
	  }
	}

	phi[5][i][j][l]*=(activity);

      }
    }
  }
  
  //clearing the memory
  destroy_4d_double_array(qB4);
  destroy_3d_double_array(qintB4);



  return Q;


};
