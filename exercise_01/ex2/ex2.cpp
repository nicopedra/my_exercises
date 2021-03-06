#include "funz.cpp"

using namespace std;

 int main () {
//####################################################################################
 Random rnd = random_initialization(1);
//####################################################################################

 int it=10000;//number of iterations
 vector<int> N = {1,2,10,100};//to decide how many random number sum in each iteration
 vector<vector<double>> hs(N.size());//matrices with 4 columns and it rows
 vector<vector<double>> he(N.size());//each column is the histogram 
 vector<vector<double>> hl(N.size());//for each value of N
 
 // hs --> standard dice
 // he --> exponential dice
 // hl --> lorentzian/cauchy dice

 map< string, vector<vector<double>> > hist = {{"hist_stand.txt",hs},{"hist_exp.txt",he},{"hist_lor.txt",hl}};
 //map that associates to each histogram the file containing its own data
 
 for (auto & el : hist) fill_hist(rnd,it,el.second,N,el.first);//function that fills histograms
 
 system("mv hist_* ../");//moving all in the right directory

 rnd.SaveSeed();//save seed

 return 0;
}



