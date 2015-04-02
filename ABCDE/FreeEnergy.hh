void FreeEnergy(double ****w, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz, double **chiMatrix){

  
  double  currentfE, oldfE, deltafE,oldfE_iter; 
  int     maxIter=500; 
  int     i,j,k,iter,chain,ii,jj; 
  double  precision=1.0e-2; 
  double  QMultiBlock,QHomo; 
  double  fEW, fEchi, fES; 
  double  epsilon, gamma;
  double  ***delphi;
  double  ****delW;
  double  ****newW;
  double  deltaW;
  double  fE_homo;
  double  box,msg;

  delW=create_4d_double_array(ChainType,Nx,Ny,Nz,"delW");
  delphi=create_3d_double_array(Nx,Ny,Nz,"delphi");
  newW=create_4d_double_array(ChainType,Nx,Ny,Nz,"newW");

  // Calculating the Homogenous Free Energy
  fE_homo=homogenousfE(chiMatrix,chi);

  //std::cout<<"Dis Copolymer Concentration:  "<<Phi_Copo_Dis<<"  Dis Homopolymer Concentration:   "<<Phi_Homo_Dis<<std::endl;
  
  msg=1.0;
  oldfE=1.0e2;
  std::ofstream outputFile("./RESULTS/fE.dat");
  do{
   
    WaveVectors(k_vector,dxyz);
    currentfE=0.0;
    deltafE=0.0;
  
    iter=0;  
    
    do{

      if(iter<500){
	epsilon=0.01;
	gamma=0.01;
      }else{
	epsilon=0.05;
	gamma=0.05;
      }
    
      iter++;
    
      fEW=0.0;
      fEchi=0.0;
      fES=0.0;

      QMultiBlock=ConcMultiBlock(phi,w,Ns,ds,k_vector,dxyz);
      QHomo=ConcHomo(phi,w,Ns,ds,k_vector,dxyz);
     
      Incomp(eta,phi,delphi);
      Phi_Copo_Ord/=(Nx*Ny*Nz);
      Phi_Homo_Ord/=(Nx*Ny*Nz);

      fEW=0.0;
      fEchi=0.0;
      
      deltaW=0.0;

    
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    for(ii=0;ii<ChainType;ii++){
	      newW[ii][i][j][k]=0.0;  
	    }
	  }
	}
      }

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){

	    for(ii=0;ii<ChainType;ii++){
	      for(jj=0;jj<ChainType;jj++){
	  
		newW[ii][i][j][k]+=chiMatrix[ii][jj]*(phi[jj][i][j][k]);
		fEchi+=phi[ii][i][j][k]*chiMatrix[ii][jj]*phi[jj][i][j][k]*dxyz[0]*dxyz[1]*dxyz[2];

	      }

	      if(ii==0){newW[ii][i][j][k]+=eta[i][j][k];}    //A
	      if(ii==1){newW[ii][i][j][k]+=eta[i][j][k];}    //C
	      if(ii==2){newW[ii][i][j][k]+=eta[i][j][k];}    //B1
	      if(ii==3){newW[ii][i][j][k]+=eta[i][j][k];}    //B2
	      if(ii==4){newW[ii][i][j][k]+=eta[i][j][k];}    //B3
	      if(ii==5){newW[ii][i][j][k]+=eta[i][j][k];}    //B4

	      fEW+=(newW[ii][i][j][k]*phi[ii][i][j][k]*dxyz[0]*dxyz[1]*dxyz[2]);
	      delW[ii][i][j][k]=newW[ii][i][j][k]-w[ii][i][j][k];
	      deltaW+=fabs(delW[ii][i][j][k]);
	    }
	 
	  }
	}
      }

      deltaW/=(Nx*Ny*Nz);
      fEchi/=(2.0*((Nx*dxyz[0])*(Ny*dxyz[1])*(Nz*dxyz[2])));
      fEW/=(((Nx*dxyz[0])*(Ny*dxyz[1])*(Nz*dxyz[2])));
    
      fES=QMultiBlock+activity*QHomo;
   
      currentfE=-fES-fEW+fEchi-fE_homo;

      deltafE=fabs(currentfE-oldfE_iter);

      //std::cout<<"Iter="<<iter<<"   fE="<<currentfE+fE_homo<<"   fE_homo="<<fE_homo<< "   delW=" << deltaW<<"   pCopo="<<Phi_Copo_Ord<<"   pHom="<<Phi_Homo_Ord<<std::endl;
      oldfE_iter=currentfE;

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){

	    for(chain=0;chain<ChainType;chain++){
	      w[chain][i][j][k]+=(gamma*delW[chain][i][j][k]-epsilon*delphi[i][j][k]);
	    }

	  }
	}
      }

    }while(deltaW>precision);//while(iter<maxIter);//


    outputFile <<currentfE<<" "<<fE_homo<<" "<<dxyz[0]*Nx<<" "<<dxyz[1]*Ny<<" "<<dxyz[2]*Nz<<std::endl;

    box=size_adjust(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);
 
   
    if(oldfE<currentfE){
      msg=0.0;
    }
    if(msg>0.5){
      oldfE=currentfE;
    }
    if(box_min==0){
      msg=0.0;
    }
    
  }while(msg>0.5);

  SaveData(phi,w,dxyz);

  Free_Energy=currentfE+fE_homo;
  Free_Energy_Homo=fE_homo;
  Lx=dxyz[0]*Nx;
  Ly=dxyz[1]*Ny;
  Lz=dxyz[2]*Nz;

  
  outputFile <<"Done"<<std::endl;
  outputFile.close();

  
  destroy_3d_double_array(delphi);
  destroy_4d_double_array(delW);
  destroy_4d_double_array(newW);
  
};
