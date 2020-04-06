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
  //double accept_reject(double,double,double);

  double retta();//generate random number distributed with 
  		 //p(x) = 2*(1-x)
  vector<double> random_direction_3D(double);//create a vector 
 					     //in a random direction in 3D space 
   };

Random random_initialization();

double error(vector<double>,vector<double>,int);//useful to evaluate the standard
						//deviation mean in the blocking method

double mean(vector<double>,int,int);

double dev_std_mean(vector<double>);

void print_vector(vector<double>, string);

void print_matrix(vector<vector<double>>, string);

void data_blocking(int,vector<double>, double, string);//passing the number of blocks N, 
                                                       //the vector conteining the N values calculated in each block,
                                                       //the reale value, and the name of the file 
						       //in which save the results

vector<double> uniform_integral(Random&,int,int,double);//passing a random generator,
						       //number of total throwns,
				           	       //number of blocks, and if necessary
						       //the lenght of the interval
double uniform_I(Random&,int,double);//returns and integral evaluated with uniform method

vector<double> importance_sampling(Random&,int,int,double);

double importace_sampling_I(Random&,int,double);//returns an integral evaluated with 
					       //importance sampling method with
				               //straight line distribution

vector<double> RW_cartesian(Random&,vector<double>&,double);//make a discrete RW step
							    //x,y,z direction (or more dimensions)
							    //modifies the passed vector	

vector<double> RW_domega(Random&,vector<double>&,double);//make a continuous RW step 
							 //generic direction in 3D space
							 //modifies the passed vector

double dist_vect(vector<double>,vector<double>);//return the square distance between two vectors
						//default: distance between a vector and the origin	

template <typename T>
vector<T> sum_vector(vector<T>,vector<T>);

template <typename T>
void print(vector<T>);

vector<double> last_data_from_datablocking(int,vector<double>,double);
//do the same thing as function data_blocking, but returns the last values for sum_prog and err_prog

void sqrt_variance_RW (Random&,int,int,int,string); //passing random generator, number of total simulations, number of blocks for data_blocking,
//number of steps for each simulation, file to save the data, int to decide if you want to
//generate a discrete RW (value 0) or a continue RW (value 1)
 
#endif 

