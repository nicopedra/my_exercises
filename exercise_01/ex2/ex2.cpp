#include "funz.cpp"

using namespace std;

 int main () {
//####################################################################################
 Random rnd = random_initialization();
//####################################################################################

 int it=10000;//number of iterations
 vector<int> N = {1,2,10,100};//to decide how many random number sum in each iteration
 vector<vector<double>> hs(N.size());//matrices with 4 columns and it rows
 vector<vector<double>> he(N.size());//each column is the histogram 
 vector<vector<double>> hl(N.size());//for each value in N
 map< string, vector<vector<double>> > hist = {{"hist_stand.txt",hs},{"hist_exp.txt",he},{"hist_lor.txt",hl}};

 for (auto & el : hist) fill_hist(rnd,it,el.second,N,el.first);
 
 system("mv hist_* ../");
 rnd.SaveSeed();

 return 0;
}



