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
  
  ReadInput.open("input.dat"); //Read input

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
  cout << "rcut: " << rcut << endl;
  ReadInput >> delta;//guardo slide per capire le unità di misura. viene usato l'epsilon dell'argon
  ReadInput >> nstep;
  ReadInput >> iprint;//ogni quanto stampare a che punto sono della simulazione

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements   //they're just indices
  //iv = 0; //Potential energy
  //ik = 1; //Kinetic energy
  //ie = 2; //Total energy
  //it = 3; //Temperature
  vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
  cout << "vtail: " << vtail << endl;
  cout << "ptail: " << ptail << endl;
  n_props = 5; //Number of observables, already add pressure
  bin_size = (box*0.5)/nbins;
  cout << "size of each bin: " << bin_size << endl;
  for (uword i=0;i<=nbins;i++) edges[i] = i*bin_size; 
  mat v(nstep/10.,nbins); //ogni riga è una gdir
  stima_gdir = v;

string start_file;
if(restart == 1) start_file = "config.0";
else start_file = "config.fcc";
//Read initial configuration
  cout << "Read initial configuration from file "+start_file << endl << endl;
  X.load(start_file);
  X = X*box;

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
else{
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


Measure(0);

//Print initial values for the potential energy and virial
  cout << "Initial potential energy (with tail corrections) = " << stima_pot+vtail << endl;
  cout << "Pressure (with tail corrections) = " << stima_press+ptail*(double)npart/vol << endl << endl;


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

void Measure(int nconf){ //Properties measurement
  double v, w, t;
  int bin;

  v = 0.0; //reset observables
  w = 0.0;
  t = 0.0;
 
 mat D_vec;
 mat Dr;
 uvec ind;
 //uvec h(nbins+1,fill::zeros);
 rowvec h(nbins*10,fill::zeros); //per sicurezza metto tante coordinate. ma quelle che mi interessano solo le prime nbins
 for (uword i=0;i<npart-1;i++) {
	 D_vec = X(span(i+1,npart-1),span::all);
	 D_vec.each_row() -= X.row(i);
	 D_vec = Pbc(-D_vec);
	 Dr = sqrt(sum(pow(D_vec,2),1));
	  
	 for (auto& el : Dr) {
		 bin = int(el/bin_size);
		 h[bin] = h[bin] + 2;
	 }
	 
	 //h += histc(Dr,edges); //in armadillo
	 ind = find(Dr<rcut && Dr>0);
	 for(auto& el : ind) {
		v += 4.0/pow(Dr[el],12) - 4.0/pow(Dr[el],6);//potential
		w += 16.0/pow(Dr[el],12) - 8.0/pow(Dr[el],6);//virial
	 }
 }
   double deltaVr;
   //vec hv =conv_to<vec>::from(h);
    for (uword i=0;i<nbins;i++){
	deltaVr = -4.*M_PI/3. * ( pow(edges(i),3) - pow(edges(i+1),3));
    	h(i) /= (rho*npart*deltaVr);
    }
    stima_gdir.row(nconf) = h(span(0,nbins-1)); 	
    t = 0.5*sum(sum(pow(V,2)));

    stima_pot = v/(double)npart;  //+vtail; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    stima_press = rho*stima_temp+(w/vol);  //+ptail*(double)npart/vol;//pressure

    //saving here to do data_blockig later
    properties[0].push_back(stima_pot+vtail);
    properties[1].push_back(stima_kin);
    properties[2].push_back(stima_temp);
    properties[3].push_back(stima_etot+vtail);
    properties[4].push_back(stima_press+ptail*(double)npart/vol);

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
  //cout << "mean temperature: " << mean(properties[2],properties[2].size())<< endl;
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

void ConfFinal(void){ //Write finals configuration

  fstream WriteOldConf("config.final",ios::out);
  cout << "Print penultimate configuration in config.final " << endl << endl;
  mat Y = Xold/box;//it is not yet the last step
  Y.print(WriteOldConf);
  WriteOldConf.close();

  fstream WriteFinalConf("config.0",ios::out);
  cout << "Print final configuration in config.0 " << endl << endl;
  X = X/box;
  X.print(WriteFinalConf);
  WriteFinalConf.close();

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
  vector<string> mark(npart,"H ");
  for (int i=0;i<npart;i++){
	  //AllXYZ << mark[i] <<"\t" << Y(i,0) <<"\t" << Y(i,1) <<"\t" << Y(i,2) << endl;
  AllXYZ <<mark[i] <<"\t"<< Y(i,0) <<"\t" << Y(i,1) <<"\t" << Y(i,2) << endl;
 
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

void data_blocking_MD(int N) {

int L = (nstep/10.)/N; //cause I measure properties each 10 steps
vector<string> names = {"ave_epot","ave_ekin","ave_temp","ave_etot","ave_press"};
string gdir_name = "output.gave.out";
int j=0;
vector<double> v_mean;
vector<double> data(2);
 for (auto & el : names) {
	for (int i=0;i<N;i++) 
		 v_mean.push_back( mean(properties[j], (i+1)*L, i*L ));
	 if ( j== 2) data = last_data_from_datablocking(N,v_mean);  
         data_blocking(N,v_mean,0,el+to_string(nstep)+".out");
	 properties[j].clear();
	 j++;
	 v_mean.clear();
 }

 accettazione = data[1]*0.8;
 m_temp = data[0];
 cout << "temperatura di ora: " << data[0] << " , con incertezza: " << data[1]<< endl;
 // radial correlation function analysis
 ofstream Gave(gdir_name,ios::out);
 vector<double> appo;
 for (int i=0;i<nbins;i++) {
 	for (auto& el : stima_gdir.col(i)) appo.push_back(el);
        for (int j=0;j<N;j++) 
		 v_mean.push_back(mean(appo, (j+1)*L, j*L ));
        data = last_data_from_datablocking(N,v_mean);
	v_mean.clear();
	appo.clear();
	Gave << (bin_size*0.5+i*bin_size) <<"\t"<<data[0] <<"\t" <<data[1]<< endl;
	v_mean.clear();
	appo.clear();
 }

Gave.close(); 
};

vector<double> last_data_from_datablocking(int N,vector<double> simulation_value) {

 vector<double> err_prog;
 vector<double> sum_prog(N,0.);
 vector<double> simulation_value2;
 vector<double> su2_prog(N,0.);
 vector<double> data(2);
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
	data = {sum_prog[N-1],err_prog[N-1]};

	return data;
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
        return sum/(double)(last_index-first_index);
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
