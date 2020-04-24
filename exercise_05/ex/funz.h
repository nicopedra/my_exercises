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

double a_0 = 0.0529e-9;

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
  vec Rannyu3D_center(vec,double);//T(x|x') uniform, centered in vec
  vec Rannyu3D_gauss(vec,double);//T(x|x') gaussian
  double Gauss(double mean, double sigma);

  double retta();//generate random number distributed with 
  		 //p(x) = 2*(1-x)
  vector<double> random_direction_3D(double);//create a vector 
 					     //in a random direction in 3D space 
  double exponential_dist(double);//generate a random number with
  				  //exponential distribution
  double lorentzian_dist(double,double);//generate random number with
  					//lorentzian distribution
   };

Random random_initialization();

double error(vector<double>,vector<double>,int);//useful to evaluate the standard
						//deviation mean in the blocking method

double mean(vector<double>,int,int);

void print_vector(vector<double>, string);

void print_matrix(vector<vector<double>>, string);

void data_blocking(int,vector<double>, double, string);//passing the number of blocks N, 
                                                       //the vector conteining the N values calculated in each block,
                                                       //the reale value, and the name of the file 
						       //in which save the results

double dist_vect(vector<double>,vector<double>);//return the square distance between two vectors
						//default: distance between a vector and the origin	
double dist(vector<double>,vector<double>);//sqrt(dist_vect)

template <typename T>
vector<T> sum_vector(vector<T>,vector<T>);

template <typename T>
void print(vector<T>);

vector<double> last_data_from_datablocking(int,vector<double>,double);
//do the same thing as function data_blocking, but returns the last values for sum_prog and err_prog

double psi_100_square (vec);//return the value of hydrogen GS distribution function

double psi_210_square (vec);//return the value of hydrogen psi_2,1,0 excited state distribution function

double A_psi100 (vec, vec);//min of[1,psi(new)/psi(old)] 

double A_psi210 (vec, vec);

void Metropolis_uniform(Random&,int,int,double);

void Metropolis_gauss(Random&,int,int,double); 
#endif 

