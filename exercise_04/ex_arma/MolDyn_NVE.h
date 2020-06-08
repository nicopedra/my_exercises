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
#include <armadillo>    //armadillo library

using namespace arma;
using namespace std; //armadillo namespace

//to remember physical scale
/*
//######## Argon ##########
  sigma = 0.34e-9;//m
  eps_kb = 120;//kelvin
  Eps = eps_kb*(1.38e-23);//Joule
  m = 39.948;//amu

  //######## Krypton ##########
  //sigma = 0.364e-9;//m
  //eps_kb = 164;//kelvin
  //Eps = eps_kb*(1.38e-23);//Joule
  //m = 83.798;//am
*/

//parameters, observables
const int m_props=5;
const int nbins = 100;
double stima_pot, stima_kin, stima_etot, stima_temp,stima_press;
//tail corrections
double vtail;
double ptail;
//mean temp to understand when stop equilibration
double m_temp;
double accettazione;

//to save properties and do blocking average
vector<vector<double>> properties(m_props);
//for g(r)
vec edges (nbins+1);
mat stima_gdir;

//configuration
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
double bin_size;
//functions

//Initialization
void Input(void);

//move particles with Verlet integrator
void Move(void);

//print actual configurations
void print_conf(void);

//overwrite config.0 and config.final
void ConfFinal(void);

void ConfXYZ(int);

//measure physical properties
void Measure(int);

//print physical istantaneous properties on file
void print_properties(); 

void Print(vector<double>,string);

//evaluate Force between particles
mat Force();

double Pbc(double);

//Pbc for matrix object
mat Pbc(mat);

double error(vector<double>,vector<double>,int);

void data_blocking(int,vector<double>,double,string);//same as the previous exercises
vector<double> last_data_from_datablocking(int,vector<double>);

double mean(vector<double>,int,int a=0);//vector, last indice, first indice
//do the blocking analysis for the MD system
void data_blocking_MD(int); 

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
