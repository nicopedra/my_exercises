#include "funz.cpp"

using namespace std;
using namespace arma;

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization();
//###########################################################

//################################# EXERCISE 1 ##################
 int N=100;//number of blocks
 int iterations=1000000;

 Metropolis_uniform(rnd,iterations,N,a_0);
 Metropolis_gauss(rnd,iterations,N,2*a_0);

 rnd.SaveSeed();

 return 0;

}



