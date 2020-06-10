#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <armadillo>

#ifndef __funz__
#define __funz__

using namespace std;
using namespace arma;

//double a_0 = 0.0529e-9;
int accepted, attempted;

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

Random random_initialization(int);

double error(vector<double>,vector<double>,int);//useful to evaluate the standard
						//deviation mean in the blocking method

double mean(vector<double>,int,int);

void print_vector(vector<double>, string);

void data_blocking(int,vector<double>, double, string);//passing the number of blocks N, 
                                                       //the vector conteining the N values calculated in each block,
                                                       //the reale value, and the name of the file 
						       //in which save the results

template <typename T>
void print(vector<T>);

vector<double> last_data_from_datablocking(int,vector<double>);
//do the same thing as function data_blocking, but returns the last values for sum_prog and err_prog

double psi_T (double,double,double);//psi_trial

double psi_T_square (double,double,double);//square modulus

double A_psi_T (double,double,double,double);//acceptance Metropolis

double V(double);//potential

double D2_psi_T(double,double,double);//second derivative of psi_trial

double H_psi_T(double,double,double);// (H | psi > )/ psi 

void Metropolis(Random&,int,int,double);

#endif 

