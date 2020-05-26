#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <armadillo>
#include <functional>
#include <complex>

#ifndef __funz__
#define __funz__

using namespace std;
using namespace arma;

int accepted, attempted;
int iprint;
int it;
//contiene le coordinate di ogni singola città
vector<vector<double>> city_coordinates;
//per crossover
vector<int> old_child_one;
vector<int> old_child_two;
vector<int> new_child_one;
vector<int> new_child_two;
vector<int> chromo_0;
double total_fitness;
vector<double> fitness;
vector<double> range_fitness;
bool crossv;

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
  double Rannyu(double min, double max);
  double Rannyu1D_center(double,double);//T(x|x') uniform, centered in vec
  double Gauss(double mean, double sigma);
   };

Random random_initialization();

class Population {

private:
        int size; //grandezza popolazione	
	int genes;//quanti geni ci sono dentro un cromosoma, 
	           //pari al numero di città in questo caso specifico
	vector<vector<int>> chromosomes;
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
  int get_gene_chromo(int j,int i) {return chromosomes[i][j];}
  void fill_initial_population();//popolazione iniziale
  void selection(Random&);
  void order_chrom();
  void crossover(Random&);
  void mutation(Random&,vector<int>&);
  void print_population(); 
  double mean_cost();
  void print_average_path(); 
  void print_best_path();
  void update();
 };

void check_function (vector<int>);

template <class InputIt1,class T>
T square_norm (InputIt1,InputIt1,T); 

template <typename T>
void shift_vector(vector<T>&,int,int); 

template <typename T>
vector<T> minus_vec(vector<T>,vector<T>);

double cost_function(vector<int>); 
void picture_path(vector<int>,string);

struct {
        bool operator()(vector<int> a, vector<int> b) const
        {   
            return cost_function(a) < cost_function(b);
        }   
    } costLess;
struct {
        bool operator()(vector<int> a, vector<int> b) const
        {   
            return cost_function(a) > cost_function(b);
        }   
    } costMax;

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

