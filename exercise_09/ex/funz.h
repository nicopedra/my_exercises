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

int accepted, attempted;
int iprint;
int it;
//contiene le coordinate di ogni singola città
vector<vector<double>> city_coordinates;
//per crossover
vector<int> new_child_one;
vector<int> new_child_two;
//per inizializzazione e check
vector<int> chromo_0;
//per selezione
double total_fitness;
vector<double> fitness;
//controllo avvenuto crossover
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
  double Gauss(double mean, double sigma);
   };

Random random_initialization(int);

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
  void selection(Random&);//algoritmo di selezione
  void order_chrom();//ordina la popolazione in base alla cost
  void crossover(Random&);//algoritmo di crossover
  void mutation(Random&,vector<int>&);//mutazione
  void print_population(); 
  double mean_cost();
  void print_average_path(); 
  void print_best_path();
  void update();//update population 
 };

//check permutation
void check_function (vector<int>);

//for cost
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
struct {
        bool operator()(vector<int> a, vector<int> b) const
        {   
            return cost_function(a) > cost_function(b);
        }   
    } costMax;

void initialize_circon(Random&,int);

void initialize_square(Random&,int);

double mean(vector<double>,int,int);

template <typename T>
void print_vector(vector<T>);

template <typename T>
void print(vector<T>);

#endif 

