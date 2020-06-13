#include "funz.cpp"

using namespace std;
using namespace arma;

//#define square

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization(2);
//###########################################################

//################################# EXERCISE 1 ##################

 //MC step per each inverse temperature
 int N=500;
 int cities = 32;
 it = 0;
 nconf = 0;
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
 int i = 0;
 double j = 100000;
 while (i<40) {
	 j /= 2.;
	 beta.push_back(1/j);
	 i++;
 }
 
 ofstream best_beta("beta_path.txt");
  for (auto& el : beta) {
	it ++;
	//MC steps
	single_MC_SA(rnd,el,N,path);
	if (it%10 == 0) {
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



