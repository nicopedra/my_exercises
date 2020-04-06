#include "funz.cpp"

using namespace std;

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization();
//###########################################################

//################################# EXERCISE 1 ##################
 int M=10000;//number of throwns
 int N=100;//number of blocks

 data_blocking(N,uniform_integral(rnd,M,N),1.,"data_uniform_I.txt");
 data_blocking(N,importance_sampling(rnd,M,N),1.,"data_importance_I.txt");

//################################ EXERCISE 2 ##################

 //M = number of simulations
 //N = number of blocks (in order to do statistics using data blocking)
 int N_step=100;//number of steps in each RW simulation

 sqrt_variance_RW (rnd,M,N,N_step,"data_RW_discr.txt",0);
 sqrt_variance_RW (rnd,M,N,N_step,"data_RW_cont.txt",1);  

 system("mv data*  ../");
 rnd.SaveSeed();

 return 0;

}



