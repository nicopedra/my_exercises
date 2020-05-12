/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"
#include <cstdlib>

//ovviamente prendo valore finale del data blocking

using namespace std;

int main(int argc, char** argv) {

if (argc!=3) {cerr << "insert temperature and restart parameter" << endl;
              return -1;}

temp = atof(argv[1]);
restart = atoi(argv[2]);

  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  print_final_values();
  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  //ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();
  ifstream Readspin("config.final");


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
if (restart==0) {
  	for (int i=0; i<nspin; ++i)
  	{
    	if(rnd.Rannyu() >= 0.5) s[i] = 1;
    	else s[i] = -1;
  	}
	}
else {
	for (int i=0;i<nspin;i++)
		Readspin >> s[i];
}

	Readspin.close();
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double r,rr;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    energy_old =Boltzmann(s[o],o);
    energy_new =Boltzmann(-s[o],o);

    if(metro==1) //Metropolis
    {
	    r = min(1.,exp(-beta*( energy_new-energy_old ) ) );
	    rr = rnd.Rannyu();
	    attempted++;
	    if(rr<=r) { s[o] = -s[o]; 
	    		accepted++;}
    }
    else //Gibbs sampling
    {	
	   p = 1. / (1+exp( beta*( energy_new-energy_old ) ) );
           attempted++;	
	   if (rnd.Rannyu()<=p){ s[o]=-s[o];
		                 accepted++;}
    }
  }
}

double Boltzmann(int sm, int ip)//valore di spin e indice dello spin
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;//è deltaE
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;
  double u2 = 0.0,chi = 0.0; 
//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);//calcolo qui u^2
     m += s[i];     //poi dopo in averages calcolo la capacità termica
// INCLUDE YOUR CODE HERE
  }
  u2 = u*u;
  chi = m*m;
  walker[iu] = u;
  walker[im] = m;
  walker[ic] = u2;
  walker[ix] = chi;
  ofstream Ene("inst_ene.out",ios::app);
  Ene << u/(double)nspin << endl;
  Ene.close();
// INCLUDE YOUR CODE HERE
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;//quante volte ho accumulato nuova informazione
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    //cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    Heat.open("output.hc.0",ios::app);
    double appo = blk_av[ic]/blk_norm/(double)nspin;
    double appo2 = nspin*blk_av[iu]*blk_av[iu]/pow(blk_norm*(double)nspin,2) ;
    stima_c = (appo-appo2)*beta*beta; //Heat Capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    Chi.open("output.chi.0",ios::app);
    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Chi
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    //cout << "----------------------------" << endl << endl;
}

void print_final_values(void) {

    const int wd=12;	
    ofstream Ene_final("output.ene.final",ios::app),Heat_final("output.hc.final",ios::app),Mag_final("output.mag.final",ios::app),Chi_final("output.chi.final",ios::app);
    Ene_final << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
    Heat_final << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
    Mag_final << setw(wd) << temp << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
    Chi_final << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
    Ene_final.close();Heat_final.close();Mag_final.close();Chi_final.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    double error;	
    if (iblk == 1) error = 0;
    else{ 
	    double p = (sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk;
	    if (p<0) error = sqrt(abs(p));
	    else error = sqrt(p);
    }
    return error;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
