#include "funz.cpp"

using namespace std;
using namespace arma;

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization(2);
//###########################################################

//################################# EXERCISE 1 ##################
 int N=100;//number of blocks
 int iterations=1000000; //number of total throwns

 //Metropolis algorithm with uniform T(x|y)
 Metropolis_uniform(rnd,iterations,N,a_0);
 //Metropolis algorithm with gaussian T(x|y)
 Metropolis_gauss(rnd,iterations,N,2*a_0);

 system("mv campionamenti* r* ../");//move in right directory

 rnd.SaveSeed();

 return 0;

}



