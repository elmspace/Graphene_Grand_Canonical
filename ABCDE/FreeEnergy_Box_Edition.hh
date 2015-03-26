double FreeEnergy_Box_Edition(double ****w_temp, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz_temp, double **chiMatrix){

  
  double  currentfE; 
  double  QMultiBlock,QHomo; 
  double  fES; 

  WaveVectors(k_vector,dxyz_temp);
  
  currentfE=0.0;
  fES=0.0;
  
  QMultiBlock=ConcMultiBlock(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  QHomo=ConcHomo(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
    
  fES=QMultiBlock+activity*QHomo;
  
  currentfE=-fES;
  return currentfE;
  
  
};
