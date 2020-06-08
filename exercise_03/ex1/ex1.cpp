#include "funz.cpp"

using namespace std;

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization(5);
//###########################################################

//################################# EXERCISE 1 ##################
 
 int n=100; //number of blocks
 int M=10000; //total generated numbers
 double S0=100.;
 double K=100.;
 double T=1.;
 double r=0.1;
 double sigma= 0.25;
 int n_step = 100; //step from 0 to T

 vector<double> v = black_scholes_analitica(S0,K,T,r,sigma);

 direct_black_scholes(rnd,M,n,S0,K,T,r,sigma,v[0],v[1]);
 discret_black_scholes(rnd,M,n,n_step,S0,K,T,r,sigma,v[0],v[1]);  

 system("mv di*  ../");
 rnd.SaveSeed();

 return 0;

}



