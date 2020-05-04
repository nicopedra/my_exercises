#include "MolDyn_NVE.h"

using namespace arma;
using namespace std;

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //seed = 1;    //Set seed for random numbers
  //srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> restart;
  if(restart == 1) 
	  cout << "reading configurations from precedent simulation: " << endl;
  else  
	  cout << "using method of random velocities: " << endl;
  ReadInput >> temp;
  cout << "target temperature = " << temp << endl;
  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;
  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;//unità sigma

  ReadInput >> rcut;
  ReadInput >> delta;//guardo slide per capire le unità di misura. viene usato l'epsilon dell'argon
  ReadInput >> nstep;
  ReadInput >> iprint;//ogni quanto stampare a che punto sono della simulazione

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements   //they're just indices
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 5; //Number of observables, already add pressure


string start_file;
if(restart == 1) start_file = "config.0";
else start_file = "config.fcc";
//Read initial configuration
  cout << "Read initial configuration from file "+start_file << endl << endl;
  X.load(start_file);
  X = X*box;

  print_conf();

if(restart == 1) {
	Xold.load("config.final");//penultimate configuration
	Xold = Xold*box;
	double sumv2_arma=0.0, fs_arma;
	mat appo = Xold;
	Move();
	V = Pbc(X-appo)/(2.0*delta);
	sumv2_arma = sum(sum(pow(V,2))); 
	sumv2_arma /= (double)npart;
	fs_arma = sqrt(3*temp/sumv2_arma);
	V = V*fs_arma;
	Xold = Pbc(X-V*delta);
     }

else {
//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   V.randu(npart,3);
   V = V-ones(npart,3)*0.5;
   mat sumv_arma = mean(V);
   for (uword i=0;i<V.n_cols;i++) V.col(i) = V.col(i)-ones(npart,1)*sumv_arma(i);
   double sumv2_arma = sum(sum(pow(V,2)));
   sumv2_arma /= (double)npart;
   double fs_arma = sqrt(3 * temp / sumv2_arma);
   V = V*fs_arma;
   Xold = Pbc(X-V*delta);
}
   return;
}


void Move(void){ //Move particles with Verlet algorithm

   mat F = Force();//parte difficile fare Force
   mat Xnew = Pbc( 2.0*X-Xold+F*(delta*delta) );
   V = Pbc(Xnew-Xold)/(2.0*delta);
   Xold = X;
   X = Xnew;

  return;
}

mat Force() {
     mat F; F.zeros(npart,3);
     mat D_vec;
     mat Dr;
     uvec ind;
     for(uword ip=0;ip<npart;ip++) {
	D_vec = X;
	D_vec.each_row() -= X.row(ip); 
	D_vec = Pbc(-D_vec);
	Dr = sqrt(sum(pow(D_vec,2),1));
	ind = find(Dr<rcut && Dr>0); 
	for(auto& el : ind) F.row(ip) = F.row(ip) + D_vec.row(el) * (48.0/pow(Dr(el),14) - 24.0/pow(Dr(el),8));
     }
     return F;
}

void Measure(){ //Properties measurement
  double v, w, t;

  v = 0.0; //reset observables
  w = 0.0;
  t = 0.0;
 
 mat D_vec;
 mat Dr;
 uvec ind;

 for (uword i=0;i<npart-1;i++) {
	 D_vec = X(span(i+1,npart-1),span::all);
	 D_vec.each_row() -= X.row(i);
	 D_vec = Pbc(-D_vec);
	 Dr = sqrt(sum(pow(D_vec,2),1));
	 ind = find(Dr<rcut && Dr>0);
	 for(auto& el : ind) {
		v += 4.0/pow(Dr[el],12) - 4.0/pow(Dr[el],6);//potential
		w += 16.0/pow(Dr[el],12) - 8.0/pow(Dr[el],6);//virial
	 }
 }	 
   t = 0.5*sum(sum(pow(V,2)));

    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_press = (rho*stima_temp+(w/vol)) / double(npart);//pressure

    //##################### RAGIONAMENTO CALCOLO PRESSIONE ###################
    //P = rho*T*kb + W/(3*V) 
    //in unità riscalate P' = rho'*T' + W'/(3*V') , dove con ' ho denotate le 
    //grandezze in unità LJ, ricordo che rho' = N/V'
    //ora penso a come calcolare il lavoro fatto dal sistema = P'V' (in modo semplice)
    //P'V' = T'*N + W'/3--> così non è ancora per unità di particella
    //mi manca dividere per N

    //saving here to do data_blockig later
    properties[0].push_back(stima_pot);
    properties[1].push_back(stima_kin);
    properties[2].push_back(stima_temp);
    properties[3].push_back(stima_etot);
    properties[4].push_back(stima_press);

    return;
}

void print_properties() {

  string name = "output_epot"+ to_string(nstep)+".dat";
  Print(properties[0],name);
  name = "output_ekin"+ to_string(nstep)+".dat";
  Print(properties[1],name);
  name = "output_temp"+ to_string(nstep)+".dat";
  Print(properties[2],name);
  name = "output_etot"+ to_string(nstep)+".dat";
  Print(properties[3],name);
  name = "output_press"+ to_string(nstep)+".dat";
  Print(properties[4],name);
  cout << "mean temperature: " << mean(properties[2],properties[2].size())<< endl;
}

void Print(vector<double> v, string name) {
   fstream fd; fd.open(name,ios::app);
   for (auto& el : v) fd << el << endl;
   fd.close();
}

void print_conf(void) {

	cout << "print actual configuration: " << endl;
	mat Y = X/box;
	Y.print();
}

void PenultimateConf(void) {

  fstream WriteConf;
  cout << "Print penultimate configuration in config.final " << endl << endl;
  WriteConf.open("config.final",ios::out);
  mat Y = X/box;//it is not yet the last step
  Y.print(WriteConf);
  WriteConf.close();
  
  return;
}

void ConfFinal(void){ //Write final configuration
  fstream WriteConf;

  cout << "Print final configuration in config.0 " << endl << endl;
  WriteConf.open("config.0",ios::out);
  X = X/box;
  X.print(WriteConf);
  WriteConf.close();

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  fstream AllXYZ;
  AllXYZ.open("traj.xyz",ios::app);
  if (nconf == 1) {
  AllXYZ << npart << endl;
  AllXYZ << "This is only a comment!" << endl;
  }
  mat Y = Pbc(X);
  vector<string> mark(npart,"LJ ");
  for (int i=0;i<npart;i++){
	  AllXYZ << mark[i] << Y(i,0) << Y(i,1) << Y(i,2) << endl;
  }
  AllXYZ.close();
}

double Pbc(double r) {  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

mat Pbc(mat M) {
 mat Y = M;
 for (auto& el : Y) el = el - box*rint(el/box);
 return Y;
}

void data_blocking_MD(int N,double sigma,double eps_kb,double Eps) {

int L = (nstep/10.)/N; //cause I measure properties each 10 steps
vector<string> names = {"ave_epot","ave_ekin","ave_temp","ave_etot","ave_press"};
int j=0;
vector<double> v_mean;

 for (auto & el : names) {
	for (int i=0;i<N;i++) 
		 v_mean.push_back( mean(properties[j], (i+1)*L, i*L ));
		 if(Eps!=0) {
       		 	if(j==2) for(auto& el : v_mean) el = el*eps_kb;//temp
			else if (j==4) for(auto& el : v_mean) el = el*Eps/(pow(sigma,3));//pressure
	 	 	else for(auto& el : v_mean) el = el*Eps;//energies 
		 }		
	 data_blocking(N,v_mean,0,el+to_string(nstep)+".out");
	 j++;
	 v_mean.clear();
 }

};

void data_blocking(int N,vector<double> simulation_value, double real_value, string file) {

 vector<double> err_prog;
 vector<double> sum_prog(N,0.);
 vector<double> simulation_value2;
 vector<double> su2_prog(N,0.);

 for (int i=0;i<N;i++) simulation_value2.push_back(simulation_value[i]*simulation_value[i]);

 for (int i=0; i<N; i++) {
         for (int j=0; j<i+1; j++) {
                 sum_prog[i] += simulation_value[j];
                 su2_prog[i] += simulation_value2[j];
         }
         sum_prog[i]/=(i+1);
         su2_prog[i]/=(i+1);
         err_prog.push_back(error(sum_prog,su2_prog,i));
 }

         fstream fd;
         fd.open(file,ios::out);
         for (int i=0; i<N;i++) fd << sum_prog[i]-real_value<<" "<< err_prog[i] << endl;
         fd.close();

};

double error(vector<double> AV, vector<double> AV2, int i) {
        if (i==0) return 0;
        else return sqrt( (AV2[i]-AV[i]*AV[i]) / double(i) );
};

double mean(vector<double> v,int last_index, int first_index) {
	double sum = 0;
	for (int i=first_index; i<last_index; i++) sum += v[i];
        return sum/(last_index-first_index);
}; 

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
