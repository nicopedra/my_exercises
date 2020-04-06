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

double Random :: Gauss(double mean, double sigma) {
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

vector<double> Random :: random_direction_3D(double step = 1.) {

        double theta = acos(1-2*Rannyu());
        double phi = Rannyu(0.,2*M_PI);
	vector<double> x;
	x.push_back(step*sin(theta)*cos(phi));
	x.push_back(step*sin(theta)*sin(phi));
	x.push_back(step*cos(theta));

 return x;
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

double mean(vector<double> v,int last_index, int first_index = 0) {
	double sum = 0;
	for (int i=first_index; i<last_index; i++) sum += v[i];
        return sum/(last_index-first_index);
};       

double dev_std_mean(vector<double> v) {
	double Mean = mean(v,v.size());
	double sum=0;
	for (auto el : v) sum+= (el-Mean)*(el-Mean);
	return sqrt(sum)/v.size();
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

vector<double> uniform_integral(Random& rnd,int M,int N,double lint=1.) {

	vector<double> integral;
	vector<double> x;
	int L = M/N;
     	for (int i=0;i<N;i++) 
		integral.push_back(uniform_I(rnd,L,lint));
return integral;//returns a vector conteining N estimates of the integral for each block
		//using the uniform method
};

double uniform_I(Random& rnd, int M,double lint=1.) {
 	double sum=0;
	for (int i=0;i<M;i++) sum+=M_PI/2. * cos(M_PI/2. *rnd.Rannyu(0.,lint));
 return sum/M*lint;
};

//devo pensare a una funzione che può assomigliare a cos(x)
//nell'intervallo che sto guardando. va bene tipo PI/2*(1-x)
//devo pensare al disegno della mia funzione!
vector<double> importance_sampling(Random& rnd,int M,int n, double lint=1.) {

	vector<double> integral;
	vector<double> x;
	int L=M/n;
     	for (int i=0;i<n;i++) {
	integral.push_back(importace_sampling_I(rnd,L,lint));
	}
return integral;
};

double importace_sampling_I(Random& rnd, int M,double lint=1.) {
 	double sum=0;
	double x;
	for (int i=0;i<M;i++){
		x = rnd.retta();
		sum+=M_PI/4.*cos(M_PI/2. *x)/(1-x);
	}
 return sum/M*lint;
};

vector<double> RW_cartesian(Random& rnd, vector<double>& x,double step = 1) {
	int coord = int(rnd.Rannyu(0,3));
	if (rnd.Rannyu() < 0.5) step = -step;
	x[coord] = x[coord]+step;
	return x;
};

vector<double> RW_domega_3D(Random& rnd, vector<double>& x,double step = 1) {

	if (rnd.Rannyu() < 0.5) step = -step;
	x = sum_vector(x,rnd.random_direction_3D(step));

return x;	

};

double dist_vect(vector<double> x,vector<double> y = {0,0,0}) {
          
	double sum=0;
        if(x.size() != y.size()) return -1;
	else for (int i=0; i<x.size();i++) sum+=(x[i]-y[i])*(x[i]-y[i]);
	return sum;
};

template <typename T>
vector<T> sum_vector(vector<T> x,vector<T> y) {

	vector<T> sum;
	if(x.size() != y.size()) cout << "error: they must have the same dimension" << "\n";
	else for (int i=0;i<x.size();i++) sum.push_back(x[i]+y[i]);

	return sum;
};

template <typename T>
void print(vector<T> v) {
	for (auto el : v) cout << el << "\n";
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

void sqrt_variance_RW (Random& rnd,int M,int N,int N_step, string file, int discr_or_cont) {
 
 vector<double> salva_dist;
 int L=M/N;
 vector<double> origin(3,0.);//start from the origin
 vector<vector<double>> v(N_step);//matrix with N_step columns, each column contains M different values of
 				  //the single step calculated for each simulation
 vector<double> v_mean(N);//vector used to do the last_data_blocking for each step
 vector<double> data;

 if(discr_or_cont==0) {
 	for (int i=0;i<M;i++) {//iterate over all simulation
        	for(int j=0;j<N_step;j++) {//doing N_step in each simulation
                	salva_dist.push_back( dist_vect( RW_cartesian(rnd,origin) ) );
			//saving square distances in a vector
        	}
         origin = {0.,0.,0.};//each simulation starts from the origin point
	 }
 }
 else {
	for (int i=0;i<M;i++) {
        	for(int j=0;j<N_step;j++) {
                	salva_dist.push_back( dist_vect( RW_domega_3D(rnd,origin) ) );
        	}
        origin = {0.,0.,0.};
	}
 }

 fstream fd;
 fd.open(file,ios::out);//open file to save the square root variance values and their uncertainties
 			//sqrt(<r^2>)
 for (int j=0;j<N_step;j++) {
        for (int i=0;i<M;i++) v[j].push_back(salva_dist[i*(N_step)+j]);//storing in a vector the
	//relative square distance of a single step for each iteration
                for (int k=0;k<N;k++) v_mean.push_back(sqrt( mean(v[j], (k+1)*L, k*L ) ));
        data = last_data_from_datablocking(N,v_mean);
        fd << data[0] << " " <<data[1] << endl;
        v_mean.clear();
        data.clear();
 }
 fd.close();
};

