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

double Random :: Rannyu1D_center(double x,double delta = 1.) {
   double x_new;
   x_new = x + Rannyu(-1.,1.)*delta;
   return x_new;
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

double Random :: exponential_dist(double lambda = 1) {  
                                                  
        double y = Rannyu();
        return -1/lambda*log(1-y);
};

double Random :: lorentzian_dist(double mu = 0 ,double gamma = 1) {
                                                          
        double y = Rannyu();
        return gamma*tan(M_PI*(y-0.5))+mu;
};

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

vector<double> last_data_from_datablocking(int N,vector<double> simulation_value) {

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

// è una somma di due gaussiane
double psi_T (double x,double mu,double sigma) {
	double appo = exp(-(x-mu)*(x-mu)/(2.0*sigma*sigma) )+exp(-(x+mu)*(x+mu)/(2.0*sigma*sigma) );
	return appo;

};

double psi_T_square (double x,double mu,double sigma) {
	return psi_T(x,mu,sigma)*psi_T(x,mu,sigma);

};

double A_psi_T (double init_x,double new_x,double mu,double sigma) {
        double appo1 = psi_T_square(new_x,mu,sigma);
	double appo2 = psi_T_square(init_x,mu,sigma);
	return min(1.,appo1/appo2);

};

//potenziale
double V(double x) {
	return x*x*x*x-2.5*x*x; 
};

//derivata seconda della psi
double D2_psi_T(double x,double mu,double sigma) {

	double e1 = exp(-(x-mu)*(x-mu)/(double)(2*sigma*sigma) ) / (sigma*sigma);
	double e2 = exp(-(x+mu)*(x+mu)/(double)(2*sigma*sigma) ) / (sigma*sigma);
	double appo1 = (1.-(x-mu)*(x-mu)/(sigma*sigma));
	double appo2 = (1.-(x+mu)*(x+mu)/(sigma*sigma));
	return - e1*appo1 - e2*appo2;

};

// h tagliato è 1
// m è 1
double H_psi_T(double x,double mu,double sigma) {
	double v = V(x);
	double h = (0.5 * D2_psi_T(x,mu,sigma) ) / (psi_T(x,mu,sigma));
	return v - h;

}

void Metropolis(Random& rnd,int iterations,int N,double delta) {

 int L = iterations/N;	
 mat campionamenti; campionamenti.zeros(iterations,1);
 vector<double> H;
 vector<double> H_mean;
 fstream fd;
 double mu0 = 1; //mi servono punti iniziali possibili
 double sigma0 = 1;
 double mu,sigma;
 double error;
 double init_x=0;// il punto è scegliere un punto intelligente per iniziare la simulazione
 double new_x;
 vector<double> data(2);
 vector<double> parameters(2);
 double H_new;
 double H_old;

 H_old = H_psi_T(init_x,mu0,sigma0);
 parameters[0] = mu0; parameters[1] = sigma0;

 //guardando la figura in python per il GS so dove guardare per il valore di mu e sigma
 //mu sicuramente tra 0.5 e 1.. sigma sicuramente non meno di 0.3 e non più di 1

 //ricerca del minimo di <H(mu,sigma)>_psitrial
 for (int l=0;l<15;l++){
 	for (int k=0;k<15;k++) {
		sigma = sigma0*( rnd.Rannyu()+0.5);//vario sigma e mu a partire da punti ragionevoli
		mu = mu0*rnd.Rannyu();
		attempted =0; accepted =0;
		init_x = 0;
		//Metropolis
 		for (uword i=0;i<iterations;i++) {
 			new_x = rnd.Rannyu1D_center(init_x,delta);
 				if (rnd.Rannyu() < A_psi_T(init_x,new_x,mu,sigma)) { 
					init_x = new_x;
					accepted++;
				}
			attempted++;
			H.push_back(H_psi_T(init_x,mu,sigma));

 		}
		//fine metropolis
		//data blocking
 		for (uword i=0;i<N;i++) 
			 H_mean.push_back( mean(H, (i+1)*L, i*L )); 
 		data = last_data_from_datablocking(N,H_mean);
		H_new = data[0];
		if(H_new < H_old) {//ricerca del minimo e salvo i parametri 
			H_old = H_new;
			error = data[1];
			parameters[0] = mu; parameters[1] = sigma;
		}
		//ripulisco per ricominciare
		H.clear();
		H_mean.clear();
	}
}

cout << endl;
cout << endl;

cout << "#########################################################" << endl;
cout << "E_min: " << H_old << endl;
cout << "con errore= " << error << endl;
ofstream Para("par.txt");
Para << parameters[0] << " " << endl << parameters[1] << endl;
cout << "mu = " << parameters[0] << " , sigma = " << parameters[1] << endl;
// ORA USO I PARAMETRI TROVATI PER ANALIZZARE MEGLIO L'ENERGIA TROVATA 
		cout << endl;
		sigma = parameters[1];
		mu = parameters[0];
		attempted =0; accepted =0;
		init_x = 0;
 		for (uword i=0;i<iterations;i++) {
 			new_x = rnd.Rannyu1D_center(init_x,delta);
 				if (rnd.Rannyu() < A_psi_T(init_x,new_x,mu,sigma)) { 
					init_x = new_x;
					accepted++;
				}
			campionamenti.row(i) = init_x;
			attempted++;
			H.push_back(H_psi_T(init_x,mu,sigma));

 		}
		cout << " rate accepted " << accepted/double(attempted) << endl;
 		for (uword i=0;i<N;i++) 
			 H_mean.push_back( mean(H, (i+1)*L, i*L )); 
 		data_blocking(N,H_mean,0,"HT.txt");
 		data = last_data_from_datablocking(N,H_mean);
		cout << endl;
		cout << "--------------------------------------------------------------" << endl;
	        cout <<"mu = "<< mu << "\t" <<"sigma = "<< sigma << endl;	
 		cout <<"H_T = "<< data[0] <<" , error_H_T = "<< data[1] << endl;
		cout << "-------------------------------------------------------------" << endl;
		cout << endl;
 
 print_vector(H,"HTall.txt");
 campionamenti.row(iterations-1) = init_x;
 campionamenti = campionamenti;
 fd.open("campionamentiT.txt",ios::out);
 campionamenti.print(fd);
 fd.close();

};

