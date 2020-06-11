#include "funz.h"


using namespace std;
using namespace arma;

Random :: Random(){};
Random :: ~Random(){};

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 ;
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  WriteSeed.close();
  return;
};

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
};

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
};

vec Random :: Rannyu3D_center(vec x,double delta = 1.) {
   vec x_new(x.size());
   for (uword i=0;i<x.size();i++) x_new[i]= x[i]+Rannyu(-1.,1.)*delta;
   return x_new;
};

vec Random :: Rannyu3D_gauss(vec x,double delta = 1.) {
   vec x_new(x.size());
   for (uword i=0;i<x.size();i++) x_new[i]= Gauss(x[i],delta);
   return x_new;
};

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
};

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
};

Random random_initialization(int lettura) {

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      for (int i=0;i<lettura;i++)
      	Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   return rnd;

};

double error(vector<double> AV, vector<double> AV2, int i) {
	if (i==0) return 0;
	else return sqrt( (AV2[i]-AV[i]*AV[i]) / double(i) );
};

double mean(vector<double> v,int last_index, int first_index = 0) {
	double sum = 0;
	for (int i=first_index; i<last_index; i++) sum += v[i];
        return sum/(last_index-first_index);
};       

void print_vector(vector<double> v, string file) {
	fstream fd;
	fd.open(file,ios::out);
	for (auto el : v) fd << el << "\n";
	fd.close();
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

template <typename T>
void print(vector<T> v) {
	for (auto el : v) cout << el << "\n";
};

double psi_100_square (vec cart_coord) {//ho messo a_0 davanti per non scalare dopo e avere tutto adimensionale
	
	return pow(a_0,-3)/M_PI*exp(-2*norm(cart_coord)/a_0);

};

double psi_210_square (vec cart_coord) {
	
	return pow(a_0,-5)/(32*M_PI)*norm(cart_coord)*norm(cart_coord)*exp(-norm(cart_coord)/a_0)*(cart_coord[2]/norm(cart_coord))*(cart_coord[2]/norm(cart_coord));

};

//acceptance probability psi100
double A_psi100 (vec init_x, vec new_x) {
        double appo1 = psi_100_square(new_x);
	double appo2 = psi_100_square(init_x);
	return min(1.,appo1/appo2);

};

//acceptance probability psi210
double A_psi210 (vec init_x, vec new_x) {
        double appo1 = psi_210_square(new_x);
	double appo2 = psi_210_square(init_x);
	return min(1.,appo1/appo2);

};


void Metropolis_uniform(Random& rnd,int iterations,int N,double delta) {

 //number of throwns per each block
 int L = iterations/N;	
 //useful for saving positions and then plotting point in 3D space
 mat campionamenti; campionamenti.zeros(iterations,3);
 //for saving istantaneous values of <r>
 vector<double> r;
 //for data blocking, <r> per each block
 vector<double> r_mean;
 //saving stuff on file
 fstream fd;
 //initial point
 vec init_x = {0.,0.,0.};
 //to store the new move
 vec new_x;

 //equilibration
 for (uword i=0;i<1000;i++) {
	new_x = rnd.Rannyu3D_center(init_x,delta);
	if (rnd.Rannyu()<=A_psi100(init_x,new_x)) init_x = new_x;
 }

 for (uword i=0;i<iterations;i++) { //cycle over all iterations
	//saving coordinates
	campionamenti.row(i) = init_x.t();
	//new move 
 	new_x = rnd.Rannyu3D_center(init_x,delta);
		//heart of Metropolis algorithm
 		if (rnd.Rannyu()<=A_psi100(init_x,new_x)) init_x = new_x;
	//saving istantaneous (and normalized) values
	r.push_back(norm(init_x/a_0));
 }

 //blocking analysis
 for (uword i=0;i<N;i++) 
	 r_mean.push_back( mean(r, (i+1)*L, i*L )); 
 data_blocking(N,r_mean,1.5,"r100.txt");
 // print istantaneous values
 print_vector(r,"r100all.txt");
 // print
 campionamenti.row(iterations-1) = init_x.t();
 campionamenti = campionamenti/a_0;
 fd.open("campionamenti100.txt",ios::out);
 campionamenti.print(fd);
 fd.close();

 //clean all 
 r.clear();
 r_mean.clear();
 campionamenti.zeros(iterations,3);
 //new initial point for psi_210
 init_x = {a_0/2.,a_0/2.,a_0/2.};

 //equilibration
 for (uword i=0;i<1000;i++) {
	new_x = rnd.Rannyu3D_center(init_x,delta);
	if (rnd.Rannyu()<=A_psi210(init_x,new_x)) init_x = new_x;
 }

 //Metropolis sampling
 for (uword i=0;i<iterations;i++) {

	campionamenti.row(i) = init_x.t();
 	new_x = rnd.Rannyu3D_center(init_x,delta);

 		if (rnd.Rannyu()<=A_psi210(init_x,new_x)) init_x = new_x;

	r.push_back(norm(init_x/a_0));
 }

 //data blocking
 for (uword i=0;i<N;i++) 
	 r_mean.push_back( mean(r, (i+1)*L, i*L )); 
 data_blocking(N,r_mean,5,"r210.txt");

 //saving results
 print_vector(r,"r210all.txt");

 campionamenti.row(iterations-1) = init_x.t();
 campionamenti = campionamenti/a_0;
 fd.open("campionamenti210.txt",ios::out);
 campionamenti.print(fd);
 fd.close();
};

//the same as above, but with gaussian T(x|y)
void Metropolis_gauss(Random& rnd,int iterations,int N,double delta) {

 int L = iterations/N;	
 mat campionamenti; campionamenti.zeros(iterations,3);
 vector<double> r;
 vector<double> r_mean;
 fstream fd;

 vec init_x ={0.,0.,0.};
 vec new_x;

 //equilibration
 for (uword i=0;i<1000;i++) {
	new_x = rnd.Rannyu3D_gauss(init_x,delta);
	if (rnd.Rannyu()<=A_psi100(init_x,new_x)) init_x = new_x;
 }

 for (uword i=0;i<iterations;i++) {
        campionamenti.row(i) = init_x.t();
        new_x = rnd.Rannyu3D_gauss(init_x,delta);
                if (rnd.Rannyu()<=A_psi100(init_x,new_x)) init_x = new_x;
        r.push_back(norm(init_x/a_0));
 }

 for (uword i=0;i<N;i++)
         r_mean.push_back( mean(r, (i+1)*L, i*L ));
 data_blocking(N,r_mean,1.5,"r100gauss.txt");

 print_vector(r,"r100allgauss.txt");

 campionamenti.row(iterations-1) = init_x.t();
 campionamenti = campionamenti/a_0;
 fd.open("campionamenti100gauss.txt",ios::out);
 campionamenti.print(fd);
 fd.close();


 r.clear();
 r_mean.clear();
 campionamenti.zeros(iterations,3);
 init_x = {a_0/2.,a_0/2.,a_0/2.};

 //equilibration
 for (uword i=0;i<1000;i++) {
	new_x = rnd.Rannyu3D_gauss(init_x,delta);
	if (rnd.Rannyu()<=A_psi210(init_x,new_x)) init_x = new_x;
 }

 for (uword i=0;i<iterations;i++) {
        campionamenti.row(i) = init_x.t();
        new_x = rnd.Rannyu3D_gauss(init_x,delta);
                if (rnd.Rannyu()<=A_psi210(init_x,new_x)) init_x = new_x;
        r.push_back(norm(init_x/a_0));
 }

 for (uword i=0;i<N;i++)
         r_mean.push_back( mean(r, (i+1)*L, i*L ));
 data_blocking(N,r_mean,5,"r210gauss.txt");

 print_vector(r,"r210allgauss.txt");

 campionamenti.row(iterations-1) = init_x.t();
 campionamenti = campionamenti/a_0;
 fd.open("campionamenti210gauss.txt",ios::out);
 campionamenti.print(fd);
 fd.close();

};
