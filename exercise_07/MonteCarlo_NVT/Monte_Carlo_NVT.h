/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__

//Random numbers
#include "random.h"
#include <cmath>

int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];
int restart;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,stima_g,err_pot,err_press,err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];

// thermodynamical state
int npart;
double beta,temp,vol,rho,box,rcut;
//############### dimension for a specific molecule #############
  double sigma;//nm
  double eps_kb;//kelvin
  double Eps;//Joule
  double m;//amu
//##############################################


// simulation
int nstep, nblk;
double delta;
//pigreco

//######## Argon ##########
  //sigma = 0.34e-9;//m
  //eps_kb = 120;//kelvin
  //Eps = eps_kb*(1.38e-23);//Joule
  //m = 39.948;//amu

  //######## Krypton ##########
  //sigma = 0.364e-9;//m
  //eps_kb = 164;//kelvin
  //Eps = eps_kb*(1.38e-23);//Joule
  //m = 83.798;//amu
//##############################

const double pi=M_PI;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
