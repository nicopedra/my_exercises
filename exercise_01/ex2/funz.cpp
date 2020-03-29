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

double Random :: exponential_dist(double lambda = 1) { 
                                                  
        double y = Rannyu();
        return -1/lambda*log(1-y);
};

double Random :: lorentzian_dist(double mu = 0 ,double gamma = 1) {//generate random number with
                                                          	   //lorentzian distribution
        double y = Rannyu();
        return gamma*tan(M_PI*(y-0.5))+mu;
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

template <typename T>
double chiquadro(vector<T> observed, double expected) {
	vector<double> chi_i; 
	for (auto el : observed) chi_i.push_back((el-expected)*(el-expected)/expected);
	return accumulate(chi_i.begin(),chi_i.end(),0.);
};

double genera_angolo_senzaPI (Random rnd) {

	double x=rnd.Rannyu();
        double y=rnd.Rannyu();
        double theta;
	while(x*x+y*y>=1){
        	x=rnd.Rannyu();
        	y=rnd.Rannyu();
	}
        theta=asin(y/sqrt(x*x+y*y));
        return theta;
};

vector<double> esperimento_Buffon(Random rnd,int n,int M, double dist_lines,double L_needle) {

 vector<int> N_hit(n,0.);
 vector<int> N_thr(n,M/n);//number of throwns in each block 
 vector<double> pi;
 double y,theta;

 //using traslational symmetry over x and 
 //periodic symmetry over y (dist_lines) 
 //and symmetry over half circle (generiting angle between -Pi/2 and Pi/2)
 for (int i=0;i<n;i++) {
         for (int j=0;j<M/n;j++) {
                y=rnd.Rannyu(0.,dist_lines);
                theta=genera_angolo_senzaPI(rnd);
                //theta=rnd.Rannyu(0.,2*M_PI);
                if(y+L_needle*sin(theta) > dist_lines || y+L_needle*sin(theta) < 0) N_hit[i]++;
                }
	 pi.push_back(2*L_needle*N_thr[i]/(N_hit[i]*dist_lines));//storing a value of pi
	 							 //for each block
        }
 return pi;
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




