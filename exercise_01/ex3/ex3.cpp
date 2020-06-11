#include "funz.cpp"

using namespace std;

 int main () {
//####################################################################################
 Random rnd = random_initialization(1); 
//####################################################################################

 int M=10000;//number of throwns
 int n=100;//number of blocks
 double L_needle=0.5,dist_lines=1.;
 //doing buffon experiment with a needle
 //it returns a vector conteining the result per-each block
 vector<double> stuzzicadente = esperimento_Buffon(rnd,n,M,dist_lines,L_needle);

 //blocking method
 data_blocking(n,stuzzicadente,M_PI,"datipi.txt");

 system("mv datipi.txt ../");

 rnd.SaveSeed();

 return 0;

 }




