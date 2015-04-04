void SaveData(double ****phi, double ****w, double *dxyz){

  int i, j ,k;

 
  std::string xyz="./MATLAB/xyz_" + std::to_string(abs(mu_homo)) + "_.dat";
  std::string ABCD="./MATLAB/ABCD_" + std::to_string(abs(mu_homo)) + "_.dat"; 
  //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++
  std::ofstream outputFile7(xyz);
  for (i=0;i<Nx;i++){
    outputFile7<<i*dxyz[0]<<" "<<i*dxyz[1]<<" "<<i*dxyz[2]<<std::endl;
  }
  outputFile7.close();  
  std::ofstream outputFile8(ABCD);
  outputFile8 <<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<0.0<<" "<<0.0<<" "<<0.0<<" "<<0.0<<std::endl;
  for (i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){//format A, C, B1+B2+B3, B4
	outputFile8<<phi[0][i][j][k]<<" "<<phi[1][i][j][k]<<" "<<phi[2][i][j][k]<<" "<<phi[3][i][j][k]<<" "<<phi[4][i][j][k]<<" "<<phi[5][i][j][k]<<std::endl;
      }
    }
  }
  outputFile8.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  /*
  std::string xyz="./MATLAB/xyz.dat";
  std::string ABCD="./MATLAB/ABCD.dat"; 
  //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++
  std::ofstream outputFile7(xyz);
  for (i=0;i<Nx;i++){
    outputFile7<<i*dxyz[0]<<" "<<i*dxyz[1]<<" "<<i*dxyz[2]<<std::endl;
  }
  outputFile7.close();  
  std::ofstream outputFile8(ABCD);
  outputFile8 <<Phi_Copo_Ord<<" "<<Phi_Homo_Ord<<" "<<0.0<<" "<<0.0<<" "<<0.0<<" "<<0.0<<std::endl;
  for (i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){//format A, C, B1+B2+B3, B4
	outputFile8<<phi[0][i][j][k]<<" "<<phi[1][i][j][k]<<" "<<phi[2][i][j][k]<<" "<<phi[3][i][j][k]<<" "<<phi[4][i][j][k]<<" "<<phi[5][i][j][k]<<std::endl;
      }
    }
  }
  outputFile8.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files


  std::string omega="./OMEGA/DURING_RUN/omega_" + std::to_string(abs(mu_homo)) + "_.dat";
  std::ofstream outputFile6(omega);
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	outputFile6 <<i<<" "<<j<<" "<<k<< " "<<w[0][i][j][k]<<" "<<w[1][i][j][k]<<" "<<w[2][i][j][k]<<" "<<w[3][i][j][k]<<" "<<w[4][i][j][k]<< " "<<w[5][i][j][k]<<std::endl;
      }
    }
  }
  outputFile6.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
};



