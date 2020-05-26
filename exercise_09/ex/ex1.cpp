#include "funz.cpp"

using namespace std;
using namespace arma;

 int main (int argc, char *argv[]) {
//###########################################################
 Random rnd = random_initialization();
//###########################################################

//################################# EXERCISE 1 ##################

//funziona! ma ancora la selezione non va bene.. arrivo a un cammino troppo lungo

 vector<int> v = {1,1,1,1,1,1,1,1,1};
 cout << square_norm(v.begin(),v.end(),0.) << endl;

 int iterations=50000;//10'000 deve andare bene
 iprint = 1000;
 int cities = 32;
 //non mi interessa sapere il lato o il raggio, li ho messi uguale a 1
 //on a circumference
 //initialize_circon(rnd,cities);
 //inside a square 
 initialize_square(rnd,cities);

 cout << "cities coordinates in a circumference: " << endl;
 for (auto el : city_coordinates) 
	 cout << el[0] << " " << el[1] << endl;
 cout << endl;
 cout << endl;
 int size_population = cities;//cities*2;
 Population P(size_population,cities,0.1,0.7);

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
 for (auto& el : chromo_0) cout << el << " ";
 cout << endl;
 //cout << cost_function(chromo_0) << endl;
 
 cout << "filling random chromosomes in the polulation.. " << endl;
 P.fill_initial_population();
 cout << endl;
 cout << "starting minimazing " << endl;
 cout << "number of iterations: " << iterations << endl;
 P.order_chrom();
 cout << "starting population: " << endl;
 P.print_population();
 P.order_chrom();

 for (it = 0;it<iterations;it++) {
	if (it % iprint == 0){ 
		cout << " ----------------------" << endl;
		cout << "iterazione numero: " << it << endl;
        	cout << " ----------------------" << endl;
	}
	//order
	 P.order_chrom();
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
cout << "population after minimazing " << endl;
cout << endl;
P.print_population();
picture_path(P.get_chromo(0),"picture_path.txt");

 rnd.SaveSeed();

 return 0;

}



