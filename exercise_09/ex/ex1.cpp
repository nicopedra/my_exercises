#include "funz.cpp"

using namespace std;
using namespace arma;

//#define square

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization(2); // con 5 e 500'000 step trovo persino un minimo più minimo per circ
//###########################################################

//################################# EXERCISE 1 ##################

 int iterations=50000;
 iprint = 1000;
 int cities = 32;

#ifdef square
 cout << cities << " cities inside a square" << endl; 
 initialize_square(rnd,cities);
#else
 cout << cities << " cities placed on a circumference" << endl;
 initialize_circon(rnd,cities);
#endif
 int size_population = cities;
 Population P(size_population,cities,0.2,0.7);

 cout << "########################" << endl;
 cout << "population size: " << P.get_size() << endl;
 cout << "genes per chromosome: "<< P.get_genes() << endl;
 cout << "probability mutation: "<< P.get_pm()<< endl;
 cout << "probability crossover: "<< P.get_pc() <<endl;
 cout << "########################" << endl;
 cout << endl;

 	for (int i=0;i<cities;i++)
		chromo_0.push_back(i);
 cout << "starting chromosome: " << endl;
 print_vector(chromo_0);
 cout << endl;
 cout << "filling random chromosomes in the polulation.. " << endl;
 P.fill_initial_population();
 cout << endl;
 cout << "starting minimazing " << endl;
 cout << "number of iterations: " << iterations << endl;
 wall_clock timer;
 timer.tic();
 P.order_chrom();

 for (it = 0;it<iterations;it++) {
	if (it % iprint == 0) { 
		cout << " ----------------------" << endl;
		cout << "iterazione numero: " << it << endl;
        	cout << " ----------------------" << endl;
	}
	 //printing fitness
	 P.print_average_path();
	 P.print_best_path();
	 //selection
	 P.selection(rnd);
	 //crossover + mutation
	 P.crossover(rnd);
	 P.update();
 }

P.order_chrom();
P.print_average_path();
P.print_best_path();

cout << "######################################" << endl;
cout << endl;
cout << "acceptance of crossover: " << accepted/double(attempted) << endl;
cout << "best cost_function: " << cost_function(P.get_chromo(0)) << endl;

double time = timer.toc();
cout << "time spent in optimizing: " << time << endl;

#ifdef square
ofstream Time("square.time",ios::out);
#else
ofstream Time("circumference.time",ios::out);
#endif

Time << time << endl;
Time.close();

picture_path(P.get_chromo(0),"picture_path.txt");

 rnd.SaveSeed();

 return 0;

}



