
#include "./global.hh"
#include "./ABCDE/parameters.hh"
#include "./ABCDE/WaveVectors.hh"
#include "./ABCDE/omega.hh"
#include "./ABCDE/omega_alphaBN.hh"
#include "./ABCDE/omega_ZnSc.hh"
#include "./ABCDE/omega_alphaBN_bilayer.hh"
#include "./ABCDE/omega_CAC.hh"
#include "./ABCDE/solvediffeq.hh"
#include "./ABCDE/Junc_Density.hh"
#include "./ABCDE/ConcMultiBlock.hh"
#include "./ABCDE/ConcHomo.hh"
#include "./ABCDE/fEhomo.hh"
#include "./ABCDE/Incomp.hh"
#include "./ABCDE/FreeEnergy_Box_Edition.hh"
#include "./ABCDE/size_adjust.hh"
#include "./ABCDE/SaveData.hh"
#include "./ABCDE/FreeEnergy.hh"

#include "./MODS/Mod1.hh"

using namespace std;


int main(){

  double ****w;
  double ***eta;
  double ****phi;
  double *chi;
  double *f;
  double ds;
  double *Ns;
  double ***k_vector;
  double *dxyz;
  double **chiMatrix;
  int i;
  
  w=create_4d_double_array(ChainType,Nx,Ny,Nz,"w");
  eta=create_3d_double_array(Nx,Ny,Nz,"eta");
  phi=create_4d_double_array(ChainType,Nx,Ny,Nz,"phi");
  chi=create_1d_double_array(15,"chi"); // There are only 15 interactions
  f=create_1d_double_array(ChainType,"f");
  Ns=create_1d_double_array(ChainType,"Ns");
  k_vector=create_3d_double_array(Nx,Ny,Nz,"k_vector");
  dxyz=create_1d_double_array(3,"dxyz");
  chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");

  long iseed;
  time_t t;
  iseed=time(&t);
  srand48(iseed);

  input_q=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);
  transformed_q=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);
  final_q=(double*)fftw_malloc(sizeof(double)*Nx*Ny*Nz);


  forward_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,input_q,transformed_q,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,i);
  inverse_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,transformed_q,final_q,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,i);
  
 
  //________________________________________________________________________________________________________________________________


  Mod1(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);
  
  //________________________________________________________________________________________________________________________________
  

  //Destroy memory allocations------------
  fftw_free(input_q);
  fftw_free(transformed_q);
  fftw_free(final_q);
  
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);

  destroy_4d_double_array(w);
  destroy_3d_double_array(eta);
  destroy_4d_double_array(phi);
  destroy_1d_double_array(chi);
  destroy_1d_double_array(Ns);
  destroy_1d_double_array(f);
  destroy_3d_double_array(k_vector);
  destroy_1d_double_array(dxyz);
  destroy_2d_double_array(chiMatrix);
  //-------------------------------------


  return 0;
}




 //omega_alphaBN(w);
  //omega_CAC(w);
  //omega_alphaBN_bilayer(w);
  //omega_ZnSc(w);

  /*
  double fE_homo;
  std::ofstream outputFile("./RESULTS/p_ave.dat");
  do{
    fE_homo=homogenousfE(chiMatrix,chi);
    std::cout<<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<mu_homo<<std::endl;
    outputFile <<Phi_Copo_Dis<<" "<<Phi_Homo_Dis<<" "<<mu_homo<<std::endl;
    mu_homo+=0.1;
  }while(Phi_Homo_Dis<0.999);
  outputFile.close();
  std::cin>>i;
  */
