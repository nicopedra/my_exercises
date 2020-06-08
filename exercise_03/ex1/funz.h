#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>

#ifndef __funz__
#define __funz__

using namespace std;

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

  double retta();//generate random number distributed with 
  		 //p(x) = 2*(1-x)
  double exponential_dist(double);//generate a random number with
  				  //exponential distribution
  double lorentzian_dist(double,double);//generate random number with
  					//lorentzian distribution
   };

Random random_initialization(int);

double error(vector<double>,vector<double>,int);//useful to evaluate the standard

void print_vector(vector<double>, string);

void print_matrix(vector<vector<double>>, string);

void data_blocking(int,vector<double>, double, string);//passing the number of blocks N, 
                                                       //the vector conteining the N values calculated in each block,
                                                       //the reale value, and the name of the file 
						       //in which save the results

vector<double> last_data_from_datablocking(int,vector<double>,double);
//do the same thing as function data_blocking, but returns the last values for sum_prog and err_prog
double N(double);//gauss cumulative function

vector<double> black_scholes_analitica(double,double,double,double,double);
//analytical solution, parameters are S0,K,T,r,sigma

void direct_black_scholes(Random&,int,int,double,double,double,double,double,double);
//the same as below, without the number of step into divide T

void discret_black_scholes(Random&,int,int,int,double,double,double,double,double,double); 
//Passing the random generator, number of total throwns, number of blocks for data_blocking
//number of step into divide T, then the parameters S0,K,T,r,sigma then the last
//2 values are the values of cut and put-option obteined from the analytical solution
#endif 

