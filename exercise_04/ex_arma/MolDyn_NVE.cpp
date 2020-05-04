/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "funz.cpp"

using namespace std;
using namespace arma;

int main(){
 wall_clock timer; //measure of esecution time
 timer.tic();
 Input();  //Inizialization
//######## Argon ##########
  sigma = 0.34e-9;//m
  eps_kb = 120;//kelvin
  Eps = eps_kb*(1.38e-23);//Joule
  m = 39.948;//amu
//######################
  //int nconf = 1;
  //print_conf();
  int N = 100; //number of blocks, data_blocking
  for(int istep=1; istep <= nstep; ++istep) {
     Move();           //Move particles with Verlet algorithm
     if (istep ==2) print_conf();
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){//every 10 steps, avoiding correlations
        Measure();     //Properties measurement
	//ConfXYZ(nconf);//Write actual configuration in XYZ format  
        //nconf += 1;
     }
     if (istep==nstep-1) PenultimateConf(); 
  }
  ConfFinal();  //Write final configuration to restart
  print_properties(); 
  data_blocking_MD(N);// normalized values
  //data_blocking_MD(N,sigma,eps_kb,Eps); //when you want values in SI unit
  double time = timer.toc();
  cout << endl;
  cout <<"################################################################" << endl;
  cout << "REMEMBER: if want to save final and penultimate configurations" << endl;
  cout <<"in file old.0 (last one) and old.final(penultimate) do command-> make copy" << endl;
  cout <<"##################################################################" << endl;
  cout << endl;
  cout << "time passed: " << time <<"s" << endl;
  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
