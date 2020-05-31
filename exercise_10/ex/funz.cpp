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

void single_MC_SA (Random& rnd,double beta, int n,vector<int>& path) {

 	ofstream Best("best_path.txt",ios::app);
	double p;
	double L_old,L_new;
	vector<int> old_path;
	double r;
	for (int i =0;i<n;i++) {
		old_path = path;
		L_old = cost_function(path);
		Mutation(rnd,path);
		L_new = cost_function(path);
		if (L_new > L_old) p = exp(-beta*(L_new-L_old)) ;
		else p = 1;
		r = rnd.Rannyu();
		if (r > p) { 
			path = old_path;
			refused++;
		}
		attempted++;
	     if(i%10 == 0) 
		     Best << it  << "\t" << cost_function(path) << endl;
	}
	Best.close();
};	

void Mutation(Random& rnd,vector<int>& v) {

//mutazione con una certa probabilità
int size = v.size();
      //diversi tipi di mutazione
      //random
      double r = rnd.Rannyu();
      if (r < 0.25)
		random_shuffle(v.begin()+1,v.end());
      //shift di n posizioni per m città contigue
      if (r >= 0.25 && r < 0.5 ) {
	        int m = int ( rnd.Rannyu(1,size-1) );//quanti shiftare, parto sempre dal secondo, non l'ultimo
		int n = int ( rnd.Rannyu(1, size - m ) );//lo shift, con un certo limite 
		shift_vector(v,m,n);
      }
      //singolo scambio di città adiacenti
      if (r >= 0.5 && r < 0.75) {
	        int m = int (rnd.Rannyu(1,size-1) );
		int n = int (rnd.Rannyu(1,size-1) );
		while (n == m) {
			n = int (rnd.Rannyu(1,size-1) );
		}
	        swap(* (v.begin()+m) ,*( v.begin()+n) );
      }
      //inversione dell'ordine delle città in un certo range
      if (r >= 0.75) {
	        int m = int (rnd.Rannyu(1,size-1));//inizio, non può essere il primo o l'ultimo
	        int n = int (rnd.Rannyu(m+1,size));//fine, può essere anche l'ultimo
	        reverse(v.begin()+m,v.begin()+n);					
      } 
		
};



//questo fa tutte le possibili permutazioni
/*
void Population :: fill_initial_population(vector<int> chromo_0) {
 
 int i = 0;
 sort(chromo_0.begin()+1,chromo_0.end());
 do {
	 if ( rnd.Rannyu() < 0.5 ) { 
		chromosomes.push_back(chromo_0);
	 	i++;
	 }
 }
 while (next_permutation(chromo_0.begin()+1,chromo_0.end()));
};
*/

void check_function (vector<int> test) {

	if (not is_permutation(path_0.begin(),path_0.end(),test.begin())) {
		cout << "non ho ricreato una permutazione delle città!" << endl;
		exit(EXIT_SUCCESS);
	}
}

double cost_function(vector<int> chr) {

	vector<double> v(chr.size());
	vector<double> appo(city_coordinates[0].size());

	for (int k = 0; k < chr.size()-1;k++) {
			for (int i=0;i<appo.size();i++) 
				appo[i]= (city_coordinates[ chr[k] ][i]-city_coordinates[ chr[k+1] ][i] );
			v[k] = square_norm(appo.begin(),appo.end(),0.); 
	}

	for (int i=0;i<appo.size();i++) 
			appo[i] = (city_coordinates[ chr[chr.size()-1] ][i]-city_coordinates[ chr[0] ][i] );
	v[chr.size()-1] = square_norm(appo.begin(),appo.end(),0.); 

        return accumulate(v.begin(),v.end(),0.);
};

template <class InputIt1,class T>
T square_norm (InputIt1 first1,InputIt1 last1,T init) {
	while (first1 != last1) {
		init += *first1* *first1;
		++first1;
	}
	return init;
}

void picture_path(vector<int> chr,string file_name) {

	ofstream BestPath(file_name,ios::out);
	for (auto& el : chr)
		BestPath << el << "\t"  << city_coordinates[el][0] << "\t" << city_coordinates[el][1] << endl;
	BestPath << chr[0] << "\t" << city_coordinates[chr[0]][0] << "\t" << city_coordinates[chr[0]][1] << endl;	
	BestPath.close();
}

template <typename T>
void shift_vector(vector<T> &v, int m,int n) {
	for (int i = m; i > 0 ; i-- ) {
		
		for (int j = 0; j < n; j++ ) {
			swap (*(v.begin()+i+j) , *(v.begin()+i+j+1) );
		}
	}
}

Random random_initialization(int lettura = 0) {

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
	for (int k = 0;k<lettura;k++)
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

//per ora raggio 1, poi si potra modificare
void initialize_circon(Random& rnd,int cities) {
        double angle;
	complex<double> z;
	for (int i =0;i<cities;i++) {
	angle = rnd.Rannyu(0.,2*M_PI);
	z = polar(1.,angle);
	city_coordinates.push_back({real(z),imag(z)});
	}
//la circonferenza si può ottenere da un cerchio unitario
//basta usare i complex con raggio 1	

}	

void initialize_square(Random& rnd,int cities) {
	double x;
	double y;
	for (int i =0;i<cities;i++) {
	x = rnd.Rannyu();
	y = rnd.Rannyu();
	city_coordinates.push_back({x,y});
	}

}


double error(vector<double> AV, vector<double> AV2, int i) {
	if (i==0) return 0;
	else return sqrt( (AV2[i]-AV[i]*AV[i]) / double(i) );
};

double mean(vector<double> v,int last_index, int first_index = 0) {
	double sum = 0;
	for (int i=first_index; i<last_index; i++) sum += v[i];
        return sum/(last_index-first_index);
};       

template <typename T>
void print_vector(vector<T> v) {

	for (auto el : v) cout << el << " ";

	cout << endl;
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
