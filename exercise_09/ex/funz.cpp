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

//questo fa tutte le possibili permutazioni
/*
void Population :: fill_initial_population(vector<int> chromo_0) {
 
 sort(chromo_0.begin()+1,chromo_0.end());
 do {
	 if ( rnd.Rannyu() < 0.5 ) { 
		chromosomes.push_back(chromo_0);
	 }
 }
 while (next_permutation(chromo_0.begin()+1,chromo_0.end()));
};
*/

void Population :: fill_initial_population() {
 attempted = 0; accepted = 0; 
 int i = 0;
 vector<int> start_chromo = chromo_0;
 while (i<size) {
	chromosomes.push_back(start_chromo); 
 	random_shuffle(start_chromo.begin()+1,start_chromo.end());
	i++;
 }

};

void Population :: selection(Random& rnd) {

for (auto& el : chromosomes) fitness.push_back(1./cost_function(el));
	total_fitness = accumulate(fitness.begin(),fitness.end(),0.);

double r;
for (;;) {
	mother = int(rnd.Rannyu(0,size));
	r = rnd.Rannyu();
	if (r < fitness[mother]/total_fitness ) 
			break;
}
for (;;) {
	father = int(rnd.Rannyu(0,size));
	r = rnd.Rannyu();
	if (r < fitness[father]/total_fitness) 
			break;
}
//libero
fitness.clear();
};

void Population :: crossover(Random& rnd) {
	
attempted++;
double r = rnd.Rannyu();
crossv = false;

if (r < p_c) {
	crossv = true;
	accepted++;

	new_child_one = chromosomes[father];
	new_child_two = chromosomes[mother];

	if (new_child_one != new_child_two) {

		taglio = genes*0.5;
		//mi tengo i geni prima del taglio e li salvo in old_child

		copy(new_child_two.begin(),
			new_child_two.end()-taglio,
			back_inserter(old_child_one));
 		copy(new_child_one.begin(),
			new_child_one.end()-taglio,
			back_inserter(old_child_two));
	
		for (auto& el : old_child_one)
			new_child_one.erase(find(new_child_one.begin(),
					new_child_one.end(),el));
		for (auto& el : old_child_two)
			new_child_two.erase(find(new_child_two.begin(),
					new_child_two.end(),el));
			
		reverse(old_child_one.begin(),old_child_one.end());
		for (auto& el : old_child_one)
			new_child_one.insert(new_child_one.begin(),el);

		reverse(old_child_two.begin(),old_child_two.end());
		for (auto& el : old_child_two)
			new_child_two.insert(new_child_two.begin(),el);
			
		//controllo che siano giuste le permutazioni
		check_function(new_child_two);
		check_function(new_child_one);
	
		old_child_one.clear();
		old_child_two.clear();

		}
	
	this->mutation(rnd,new_child_one);
	this->mutation(rnd,new_child_two);

	}
}

void Population :: mutation(Random& rnd,vector<int>& v) {

//mutazione con una certa probabilità
double r = rnd.Rannyu();
if (r < p_m) {
      //diversi tipi di mutazione
      //random
      r = rnd.Rannyu();
      if (r < 0.25)
		random_shuffle(v.begin()+1,v.end());
      //shift di n posizioni per m città contigue
      if (r >= 0.25 && r < 0.5 ) {
	        int m = int ( rnd.Rannyu(1,genes-1) );//quanti shiftare, parto sempre dal secondo, non l'ultimo
		int n = int ( rnd.Rannyu(1, genes - m ) );//lo shift, con un certo limite 
		shift_vector(v,m,n);
      }
      //singolo scambio di città adiacenti
      if (r >= 0.5 && r < 0.75) {
	        int m = int (rnd.Rannyu(1,genes-1) );
	        swap(* (v.begin()+m) ,*( v.begin()+m+1) );
      }
      //inversione dell'ordine delle città in un certo range
      if (r >= 0.75) {
	        int m = int (rnd.Rannyu(1,genes-1));//inizio, non può essere il primo o l'ultimo
	        int n = int (rnd.Rannyu(m+1,genes));//fine, può essere anche l'ultimo
	        reverse(v.begin()+m,v.begin()+n);					
      } 
}      
		
};

void Population :: print_population() {

	for (auto& el : chromosomes) {
		print_vector(el);
		cout <<"relative cost_function values: " << cost_function(el) << endl;
		cout << endl;
	}
}

double Population :: mean_cost() {

	vector<double> cost;
	for (int i = 0;i<size*0.5;i++)
		cost.push_back(cost_function(chromosomes[i]));
	return mean(cost,cost.size(),0);
}

void Population :: print_average_path() {

	ofstream Path("path_average.txt",ios::app);
	Path << it << "\t" << this->mean_cost() << endl;
        Path.close();	

}

void Population :: print_best_path() {

	ofstream Path("best_path.txt",ios::app);
	Path << it << "\t" << cost_function(chromosomes[0]) << endl;
        Path.close();	
}

void Population :: update() {

if (crossv) {

	chromosomes.push_back(new_child_one);
	chromosomes.push_back(new_child_two);
	this->order_chrom();

	for (int i = 0;i<2;i++) {
        chromosomes.erase(chromosomes.end());
	}
    }

	new_child_one.clear();
        new_child_two.clear();

return;

}

void Population :: order_chrom() {
        sort(chromosomes.begin(), chromosomes.end(), costLess);	
}

void check_function (vector<int> test) {

	if (not is_permutation(chromo_0.begin(),chromo_0.end(),test.begin())) {
		cout << "non ho ricreato una permutazione delle città!" << endl;
		exit(EXIT_SUCCESS);
	}

}

double cost_function(vector<int> chr) {

	vector<double> v;
	vector<double> appo;

	for (int k = 0; k < chr.size();k++) {
		if (k!=chr.size()-1){ 
			appo = minus_vec(city_coordinates[chr[k]],city_coordinates[chr[k+1]]);
			//v.push_back(sqrt(square_norm(appo.begin(),appo.end(),0.))); //L1
			v.push_back(square_norm(appo.begin(),appo.end(),0.)); //L2
		}
		else {
			appo = minus_vec(city_coordinates[chr[k]],city_coordinates[chr[0]]);
			//v.push_back( sqrt(  square_norm(appo.begin(),appo.end(),0.) ) ); //L1
			v.push_back(square_norm(appo.begin(),appo.end(),0.) ) ; //L2
		}
		//cout << v[k] << endl;
	}
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

template <typename T>
vector<T> minus_vec(vector<T> w,vector<T> z) {
	if (w.size() != z.size()) {
		cout << "non si può fare" << endl;
		exit(EXIT_SUCCESS);
	}
	else { 
		int size = w.size();
		vector<T> v;
		for (int i=0;i<size;i++) 
			v.push_back(w[i]-z[i]);
		return v;
	}
}

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
