/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <armadillo>

using namespace arma;
using namespace std;

//parameters, observables
const int m_props=5;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp,stima_press;
vector<vector<double>> properties(m_props);
//############### dimension for a specific molecule #############
  double sigma;//nm
  double eps_kb;//kelvin
  double Eps;//Joule
  double m;//amu
//##############################################

// averages
double acc,att;

//configuration
const int m_part=108;
mat X;
mat Xold;
mat V;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;
int restart;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(void);
void Move(void);
void print_conf(void);
void PenultimateConf(void); 
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void print_properties(); 
void Print(vector<double>,string);
mat Force();
double Pbc(double);
mat Pbc(mat);
double error(vector<double>,vector<double>,int);
void data_blocking(int,vector<double>,double,string);//number of blocks, vector containing the means, real value, file_name to save all
double mean(vector<double>,int,int a=0);//vector, last indice, first indice
void data_blocking_MD(int,double=0,double=0,double=0); 
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
