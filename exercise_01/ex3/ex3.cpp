#include "funz.cpp"

using namespace std;

 int main () {
//####################################################################################
 Random rnd = random_initialization(); 
//####################################################################################

 int M=1000000;//number of throwns
 int n=100;//number of blocks
 double L_needle=0.5,dist_lines=1.;
 data_blocking(n,esperimento_Buffon(rnd,n,M,dist_lines,L_needle),M_PI,"datipi.txt");

 system("mv datipi.txt ../");

 rnd.SaveSeed();

 return 0;

 }




