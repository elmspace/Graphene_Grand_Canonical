void Junc_Density(double ****qA,double ****qC,double ****qB1,double ****qB2,double ****qB3,double ****qdagA,double ****qdagC,double ****qdagB1,double ****qdagB2,double ****qdagB3,double *Ns,double ds, double Q){

  int i,j,l,k;
  double      ***JB1C, ***JCB2, ***JB2A, ***JAB3;
  
  JB1C=create_3d_double_array(Nx,Ny,Nz,"JB1C");
  JCB2=create_3d_double_array(Nx,Ny,Nz,"JCB2");
  JB2A=create_3d_double_array(Nx,Ny,Nz,"JB2A");
  JAB3=create_3d_double_array(Nx,Ny,Nz,"JAB3");
  

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	JB1C[i][j][l]=qB1[i][j][l][(int)Ns[2]]*qdagB1[i][j][l][0]*ds;
	JCB2[i][j][l]=qC[i][j][l][(int)Ns[1]]*qdagC[i][j][l][0]*ds;
	JB2A[i][j][l]=qB2[i][j][l][(int)Ns[3]]*qdagB2[i][j][l][0]*ds;	
	JAB3[i][j][l]=qA[i][j][l][(int)Ns[0]]*qdagA[i][j][l][0]*ds;

	JB1C[i][j][l]*=(1.0/Q);
	JCB2[i][j][l]*=(1.0/Q);
	JB2A[i][j][l]*=(1.0/Q);
	JAB3[i][j][l]*=(1.0/Q);
      }
    }
  }

  

  //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++  
  std::ofstream outputFile38("./ABCD_Junc.dat");
  for (i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){//format A, C, B1+B2+B3, B4
	outputFile38<<JB1C[i][j][k]<<" "<<JCB2[i][j][k]<<" "<<JB2A[i][j][k]<<" "<<JAB3[i][j][k]<<std::endl;
      }}}
  outputFile38.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  


   destroy_3d_double_array(JB1C);
   destroy_3d_double_array(JCB2);
   destroy_3d_double_array(JB2A);
   destroy_3d_double_array(JAB3);

};
