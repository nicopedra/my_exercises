#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <complex>
#include "mpi.h"

#ifndef __funz__
#define __funz__

using namespace std;

int iterations;
vector<vector<double>> city_coordinates;
int N_migr;
vector<int> simple_chromo;
bool if_square;

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  //constructors
   Random();
  // destructor
   ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max=0);
  double Rannyu1D_center(double,double);//T(x|x') uniform, centered in vec
  double Gauss(double mean, double sigma);
   };

Random random_initialization(int);

class Population {

private:
	Random rnd;
	int accepted, attempted;//per controllare il crossover
	int it; //conta a che passo sono, utile per stampare
        int size; //grandezza popolazione	
	int genes;//quanti geni ci sono dentro un cromosoma, 
	           //pari al numero di città in questo caso specifico
	vector<vector<int>> chromosomes;
	vector<int> new_child_one;
	vector<int> new_child_two;
	double total_fitness;
	vector<double> fitness;
	bool crossv;
	vector<int> start_chromo;
	vector<vector<int>> the_BigFour;
	vector<double> Path;
	vector<double> Best_Path;
	double p_m;//probability mutation
	double p_c;//probability crossover 
	int mother,father; //for crossover and selection
	int taglio; //per il crossover
protected:

public:
  //constructors
   Population(int size,int genes,double p_m,double p_c) : size(size),genes(genes),p_m(p_m),p_c(p_c) {}
  // destructor
   ~Population() {}
  // methods
  double get_genes() {return genes;}
  double get_pc () {return p_c;}
  double get_pm () {return p_m;}
  void set_p_c (double prob_cross) {p_c = prob_cross;}
  void set_p_m (double prob_mut) {p_m=prob_mut;}
  int get_size() {return size;}
  void set_size(int n) {size = n;}
  vector<int> get_chromo(int i) {return chromosomes[i];}
  vector<int> get_theRealOne(int i) {return the_BigFour[i];}
  int get_gene_chromo(int j,int i) {return chromosomes[i][j];}
  void fill_initial_population(int);//popolazione iniziale
  void selection();
  void order_chrom();
  void crossover();
  void mutation(vector<int>&);
  void print_population(); 
  double mean_cost();
  void print_average_path(); 
  void print_best_path();
  void update();
  void check_function (vector<int>);
  void scambio(vector<int>&);
  void Exit();
  void best_list(vector<int>);
  void print_result(int);
 };

template <class InputIt1,class T>
T square_norm (InputIt1,InputIt1,T); 

template <typename T>
void shift_vector(vector<T>&,int,int); 

double cost_function(vector<int>); 
void picture_path(vector<int>,string);

struct {
        bool operator()(vector<int> a, vector<int> b) const
        {   
            return cost_function(a) < cost_function(b);
        }   
    } costLess;

void initialize_circon(Random&,int);

void initialize_square(Random&,int);

double error(vector<double>,vector<double>,int);//useful to evaluate the standard
						//deviation mean in the blocking method

double mean(vector<double>,int,int);

template <typename T>
void print_vector(vector<T>);

void data_blocking(int,vector<double>, double, string);//passing the number of blocks N, 
                                                       //the vector conteining the N values calculated in each block,
                                                       //the reale value, and the name of the file 
						       //in which save the results

template <typename T>
vector<T> sum_vector(vector<T>,vector<T>);

template <typename T>
void print(vector<T>);

vector<double> last_data_from_datablocking(int,vector<double>);
//do the same thing as function data_blocking, but returns the last values for sum_prog and err_prog

#endif 

