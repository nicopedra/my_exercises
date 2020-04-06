#include "funz.h"


using namespace std;

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

double Random :: Gauss(double mean=0, double sigma=1) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
};

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
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

double Random :: exponential_dist(double lambda = 1) {  
                                                  
        double y = Rannyu();
        return -1/lambda*log(1-y);
};

double Random :: lorentzian_dist(double mu = 0 ,double gamma = 1) {
                                                          
        double y = Rannyu();
        return gamma*tan(M_PI*(y-0.5))+mu;
};
/*
double Random :: accept_reject(double xmin=0, double xmax=1, double pmax=1) {

	double x = Rannyu(xmin,xmax);
	double r = Rannyu(); 
	if ( r < (p(x)/pmax ) return x;
	else return accept_reject();
};	
*/
double Random :: retta() {
	return 1.-sqrt(1.-Rannyu());
};

Random random_initialization() {

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
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

void print_vector(vector<double> v, string file) {
	fstream fd;
	fd.open(file,ios::out);
	for (auto el : v) fd << el << "\n";
	fd.close();
};

void print_matrix(vector<vector<double>> m, string file) {
	fstream fd;
	fd.open(file,ios::out);
	int n_row = m[0].size();
       	int n_col = m.size();
	for (int j=0;j<n_row;j++){
		for (int i=0; i<n_col; i++) {
			fd << m[i][j] << " ";
		}
		fd << endl;
	}
				
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

vector<double> last_data_from_datablocking(int N,vector<double> simulation_value, double real_value= 0.) {

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
  
 vector<double> data = {sum_prog[N-1],err_prog[N-1]};

	return data;
};

double N(double x) { 
	return 0.5*(1+erf( x/sqrt(2) ) );
};

vector<double> black_scholes_analitica(double S0,double K,double T,double r,double sigma) {
	double d1 = 1./(sigma*sqrt(T)) * (log(S0/K) + (r+(sigma*sigma)/2) * T);
        double d2 = d1-sigma*sqrt(T);
        double C = S0*N(d1)-K*exp(-r*T)*N(d2);
        double P = S0 * ( N(d1)-1. ) - K* exp(-r*T)*(N(d2)-1);

return {C,P};
};

void direct_black_scholes(Random rnd,int M,int n,double S0,double K,double T,double r,double sigma,double realc,double realp) {

	vector<vector<double>> CP(2);
	double appo;
	double sumc, sump;
	int L = M/n;
	for (int i=0;i<n;i++) {
		sumc=0;sump=0;
		for (int j=0;j<L;j++) { 
			appo = S0*exp((r-sigma*sigma/2.)*T + sigma*rnd.Gauss()*sqrt(T));
			sumc+=exp(-r*T)*max(0.,appo-K);
			sump+=exp(-r*T)*max(0.,K-appo);
		}
		CP[0].push_back(sumc/L);
                CP[1].push_back(sump/L);
       }
	data_blocking(n,CP[0],realc,"directC.txt");
	data_blocking(n,CP[1],realp,"directP.txt");		
};

void discret_black_scholes(Random rnd,int M,int n,int n_step,double S0,double K,double T,double r,double sigma,double realc,double realp) {

        vector<vector<double>> CP(2);
        double appo;
	int L = M/n;
	double sumc,sump;
        for (int j=0;j<n;j++) {
		sumc=0;sump=0;
                for (int i=0;i<L;i++) {
			appo = S0;
			for (int k=0;k<n_step;k++) 
                		 appo*=exp( (r-sigma*sigma/2.)*(T/n_step) + sigma*rnd.Gauss()*sqrt(T/n_step) );
			sumc+=exp(-r*T)*max(0.,appo-K);
			sump+=exp(-r*T)*max(0.,K-appo);
                }
		CP[0].push_back(sumc/L);
                CP[1].push_back(sump/L);
        }
	data_blocking(n,CP[0],realc,"discreteC.txt");
	data_blocking(n,CP[1],realp,"discreteP.txt");

};

