double ConcMultiBlock(double ****phi,double ****w,double *Ns,double ds,double ***k_vector,double *dxyz){

  int         i,j,l,s;
  double      Q;

  double      ****qA;
  double      ****qC;
  double      ****qB1,****qB2,****qB3;

  double      ****qdagA;
  double      ****qdagC;
  double      ****qdagB1,****qdagB2,****qdagB3;

  double      ***qintA;
  double      ***qintC;
  double      ***qintB1,***qintB2,***qintB3;

 
  qA=create_4d_double_array(Nx,Ny,Nz,((int)Ns[0]+1),"qA");
  qC=create_4d_double_array(Nx,Ny,Nz,((int)Ns[1]+1),"qC");
  qB1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[2]+1),"qB1");
  qB2=create_4d_double_array(Nx,Ny,Nz,((int)Ns[3]+1),"qB2");
  qB3=create_4d_double_array(Nx,Ny,Nz,((int)Ns[4]+1),"qB3");

  qdagA=create_4d_double_array(Nx,Ny,Nz,((int)Ns[0]+1),"qdagA");
  qdagC=create_4d_double_array(Nx,Ny,Nz,((int)Ns[1]+1),"qdagC");
  qdagB1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[2]+1),"qdagB1");
  qdagB2=create_4d_double_array(Nx,Ny,Nz,((int)Ns[3]+1),"qdagB2");
  qdagB3=create_4d_double_array(Nx,Ny,Nz,((int)Ns[4]+1),"qdagB3");

  qintA=create_3d_double_array(Nx,Ny,Nz,"qintA");
  qintC=create_3d_double_array(Nx,Ny,Nz,"qintC");
  qintB1=create_3d_double_array(Nx,Ny,Nz,"qintB1");
  qintB2=create_3d_double_array(Nx,Ny,Nz,"qintB2");
  qintB3=create_3d_double_array(Nx,Ny,Nz,"qintB3");

  //+++++++++++++++++++++++++++++++++++++++Forward++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB1[i][j][l]=1.0;
      }
    }
  }
  solveModDiffEqn_FFT(qB1,w[2],qintB1,ds,(int)Ns[2],1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintA[i][j][l]=qB1[i][j][l][(int)Ns[2]];
      }
    }
  }
  solveModDiffEqn_FFT(qA,w[0],qintA,ds,(int)Ns[0],1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB2[i][j][l]=qA[i][j][l][(int)Ns[0]];
      }
    }
  }
  solveModDiffEqn_FFT(qB2,w[3],qintB2,ds,(int)Ns[3],1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintC[i][j][l]=qB2[i][j][l][(int)Ns[3]];
      }
    }
  }
  solveModDiffEqn_FFT(qC,w[1],qintC,ds,(int)Ns[1],1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB3[i][j][l]=qC[i][j][l][(int)Ns[1]];
      }
    }
  }
  solveModDiffEqn_FFT(qB3,w[4],qintB3,ds,(int)Ns[4],1,k_vector,dxyz);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //+++++++++++++++++++++++++++++++++++++++++++Backward+++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB3[i][j][l]=1.0;
      }
    }
  }
  solveModDiffEqn_FFT(qdagB3,w[4],qintB3,ds,(int)Ns[4],-1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintC[i][j][l]=qdagB3[i][j][l][(int)Ns[4]];
      }
    }
  }
  solveModDiffEqn_FFT(qdagC,w[1],qintC,ds,(int)Ns[1],-1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB2[i][j][l]=qdagC[i][j][l][(int)Ns[1]];
      }
    }
  }
  solveModDiffEqn_FFT(qdagB2,w[3],qintB2,ds,(int)Ns[3],-1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintA[i][j][l]=qdagB2[i][j][l][(int)Ns[3]];
      }
    }
  }
  solveModDiffEqn_FFT(qdagA,w[0],qintA,ds,(int)Ns[0],-1,k_vector,dxyz);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB1[i][j][l]=qdagA[i][j][l][(int)Ns[0]];
      }
    }
  }
  solveModDiffEqn_FFT(qdagB1,w[2],qintB1,ds,(int)Ns[2],-1,k_vector,dxyz);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 
  //++++++++++++++++++++++++++++++++++++++Single Chain Partition Function+++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Q=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	Q+=(qdagB1[i][j][l][(int)Ns[2]])*dxyz[0]*dxyz[1]*dxyz[2];
      }
    }
  }
  // Normalizing with respect to the volume of the box
  Q/=((dxyz[0]*Nx)*(dxyz[1]*Ny)*(dxyz[2]*Nz));

  
  // Here we do the concentration calculation
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){


	phi[0][i][j][l]=0.0;  //A
	phi[1][i][j][l]=0.0;  //C
	phi[2][i][j][l]=0.0;  //B1
	phi[3][i][j][l]=0.0;  //B2
	phi[4][i][j][l]=0.0;  //B3


	//A1
	for(s=0;s<(Ns[0]+1);s++){
	  if(s==0 || s==(int)Ns[0]){
	    phi[0][i][j][l]+=0.5*qA[i][j][l][s]*qdagA[i][j][l][(int)Ns[0]-s]*ds;
	  }else{
	    phi[0][i][j][l]+=qA[i][j][l][s]*qdagA[i][j][l][(int)Ns[0]-s]*ds;
	  }
	}

	//C
	for(s=0;s<(Ns[1]+1);s++){
	  if(s==0 || s==(int)Ns[1]){
	    phi[1][i][j][l]+=0.5*qC[i][j][l][s]*qdagC[i][j][l][(int)Ns[1]-s]*ds;
	  }else{
	    phi[1][i][j][l]+=qC[i][j][l][s]*qdagC[i][j][l][(int)Ns[1]-s]*ds;
	  }
	}


	//B1
	for(s=0;s<(Ns[2]+1);s++){
	  if(s==0 || s==(int)Ns[2]){
	    phi[2][i][j][l]+=0.5*qB1[i][j][l][s]*qdagB1[i][j][l][(int)Ns[2]-s]*ds;
	  }else{
	    phi[2][i][j][l]+=qB1[i][j][l][s]*qdagB1[i][j][l][(int)Ns[2]-s]*ds;
	  }
	}

	//B2
	for(s=0;s<(Ns[3]+1);s++){
	  if(s==0 || s==(int)Ns[3]){
	    phi[3][i][j][l]+=0.5*qB2[i][j][l][s]*qdagB2[i][j][l][(int)Ns[3]-s]*ds;
	  }else{
	    phi[3][i][j][l]+=qB2[i][j][l][s]*qdagB2[i][j][l][(int)Ns[3]-s]*ds;
	  }
	}

	//B3
	for(s=0;s<(Ns[4]+1);s++){
	  if(s==0 || s==(int)Ns[4]){
	    phi[4][i][j][l]+=0.5*qB3[i][j][l][s]*qdagB3[i][j][l][(int)Ns[4]-s]*ds;
	  }else{
	    phi[4][i][j][l]+=qB3[i][j][l][s]*qdagB3[i][j][l][(int)Ns[4]-s]*ds;
	  }
	}

      }
    }
  }
  
  //Junc_Density(qA,qC,qB1,qB2,qB3,qdagA,qdagC,qdagB1,qdagB2,qdagB3,Ns,ds,Q);
  

  //clearing the memory
  destroy_4d_double_array(qA);
  destroy_4d_double_array(qC);
  destroy_4d_double_array(qB1);
  destroy_4d_double_array(qB2);
  destroy_4d_double_array(qB3);

  destroy_4d_double_array(qdagA);
  destroy_4d_double_array(qdagC);
  destroy_4d_double_array(qdagB1);  
  destroy_4d_double_array(qdagB2);  
  destroy_4d_double_array(qdagB3); 

  destroy_3d_double_array(qintA);
  destroy_3d_double_array(qintC);
  destroy_3d_double_array(qintB1);
  destroy_3d_double_array(qintB2);
  destroy_3d_double_array(qintB3);



  return Q;


};
