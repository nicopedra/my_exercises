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
   };

Random random_initialization();

double error(vector<double>,vector<double>,int);//useful to evaluate the standard
						//deviation mean in the blocking method

double mean(vector<double>,int,int);

double dev_std_mean(vector<double>);

void print_vector(vector<double>, string);

template <typename T>
double chiquadro(vector<T> observed, double expected);

void data_blocking(int,vector<double>, double, string);//passing the number of blocks N, 
						       //the vector conteining the N values calculated in each block,
						       //the reale value, and the name of the file in which save the results


#endif 

