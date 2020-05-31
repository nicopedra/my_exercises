#include "funz.cpp"

using namespace std;
using namespace arma;

#define square

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization(2);
//###########################################################

//################################# EXERCISE 1 ##################

 //MC step per each inverse temperature
#ifdef square
 int N = 5000;
#else
 int N=1000;
#endif

 int cities = 32;
 it = 0;
#ifdef square
 cout << cities << " cities inside a square " << endl; 
 initialize_square(rnd,cities);
#else
 cout << cities << " cities placed ona circumference " << endl;
 initialize_circon(rnd,cities);
#endif

 for (int i=0;i<cities;i++)
	path_0.push_back(i);
 cout << "starting path: " << endl;
 print_vector(path_0);
 cout << endl;
 cout << "starting minimazing " << endl;
 cout << "number of steps for each MC simulations: " << N << endl;
 vector<int> path = path_0;
 vector<double> beta;

 //filling inverse temperature 
 for (double i = 2.;i>=0.0002;i-=0.0002) beta.push_back(1./i);//per quadrato
 //cout << beta.size() << endl;
 //print_vector(beta);
 
 ofstream best_beta("beta_path.txt");
  for (auto& el : beta) {
	it ++;
	single_MC_SA(rnd,el,N,path);
	if (it%1000 == 0) {
		cout << "--------------------------------------" << endl;
		cout << "acceptance MC with inverse temperature " << el <<" is :" << (attempted-refused)/double(attempted) << endl;
		cout << "--------------------------------------" << endl;
	}
	best_beta << el << "\t" << cost_function(path) << endl;
	refused = 0; attempted = 0;
 }
cout << "best cost function: " << cost_function(path) << endl;
picture_path(path,"picture_path.txt");

rnd.SaveSeed();

 return 0;

}



