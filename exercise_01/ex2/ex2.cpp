#include "funz.cpp"

using namespace std;

 int main () {
//####################################################################################
 Random rnd = random_initialization();
//####################################################################################

 int it=10000;//number of iterations
 vector<int> N = {1,2,10,100};//to decide how many random number sum in each iteration
 vector<vector<double>> hist_stand(N.size());//matrices with 4 columns and it rows
 vector<vector<double>> hist_exp(N.size());//each column is the histogram 
 vector<vector<double>> hist_lor(N.size());//for each value in N
 double sum_stand,sum_exp,sum_lor;
 int l=0;
 fstream fd;

 for(auto & el : N) {//for each N fill the histograms
	for (int i=0;i<it;i++) {
		sum_stand=0;
		for (int j=0; j<el;j++) sum_stand+=rnd.Rannyu();
		hist_stand[l].push_back(sum_stand/el);//standard dice
	}
	for (int i=0;i<it;i++) {
		sum_exp=0;
		for (int j=0; j<el;j++) sum_exp+=rnd.exponential_dist();
		hist_exp[l].push_back(sum_exp/el);//exponential dice
	}
	for (int i=0;i<it;i++) {
		sum_lor=0;
		for (int j=0; j<el;j++) sum_lor+=rnd.lorentzian_dist();
		hist_lor[l].push_back(sum_lor/el);//lorentzian dice
	}
 l++;
 }


 print_matrix(hist_stand,"hist_stand.txt");
 print_matrix(hist_exp,"hist_exp.txt");
 print_matrix(hist_lor,"hist_lor.txt");

  system("mv hist_* ../");
  rnd.SaveSeed();

 return 0;
}



