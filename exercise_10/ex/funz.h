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

int refused, attempted;
//contiene le coordinate di ogni singola città
vector<vector<double>> city_coordinates;
//per crossover
vector<int> path_0;
int it;
int nconf;

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

void single_MC_SA (Random&,double,int,vector<int>&); 

void Mutation(Random&,vector<int>&); 

void check_function (vector<int>);

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

