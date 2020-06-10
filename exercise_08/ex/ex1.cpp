#include "funz.cpp"

//h tagliato = 1
//m = 1

using namespace std;
using namespace arma;

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization(1);
//###########################################################

//################################# EXERCISE 1 ##################
 int N=100;//number of blocks
 int iterations=1000000;

 Metropolis(rnd,iterations,N,2.5); //per scegliere il delta della T(x|x') ho guardato l'accettazione

 rnd.SaveSeed();

 return 0;

}



