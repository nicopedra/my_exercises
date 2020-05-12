/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "funz.cpp"

//#define equilibration

using namespace std;
using namespace arma;

int main(int argc, char** argv){

 int tentativo = 1;
 temp = 0.8;
 m_temp = 0;
 int nconf;
 int N;
 accettazione = 0.001;
 double errore = abs(m_temp - temp);
#ifdef equilibration
 cout << "equilibration phase" << endl;
 while ( errore > accettazione ) { //utile per equilibrazione
	cout <<"################################################################" << endl;
        cout << "                tentativo numero: " << tentativo << endl;
	cout <<"################################################################" << endl;
	if (tentativo == 1)
 		restart = 0;
	else
		restart = 1;

#else
if (argc!=2) {cerr << "insert restart parameter" << endl;
              return -1;}
 
 restart = atoi(argv[1]);

#endif
	cout << "restart: " << restart << endl;
 	wall_clock timer; //measure of esecution time
 	timer.tic();
 	Input();  //Inizialization
  	nconf = 1;
  	N = 100; //number of blocks, data_blocking
  	for(int istep=1; istep <= nstep; ++istep) {
     		Move();           //Move particles with Verlet algorithm
     		if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     		if(istep%10 == 0){//every 10 steps, avoiding correlations
        		Measure(nconf-1);     //Properties measurement
			//ConfXYZ(nconf);//Write actual configuration in XYZ format  
        		nconf += 1;
     		}
  	}
  	print_properties();
  	data_blocking_MD(N);// normalized values
  	double time = timer.toc();
	cout << "time passed: " << time <<"s" << endl;
	cout << endl;
	errore = abs(m_temp - temp);
	cout << "ora l'errore tra la temperatura del sistema e quella target è: " << errore << endl;
	system ("make copy");
	ConfFinal();
#ifdef equilibration
	tentativo++;
  }
#endif
  cout << endl;
  cout <<"################################################################" << endl;
  cout << "REMEMBER: if want to save final and penultimate configurations" << endl;
  cout <<"in file old.0 (last one) and old.final(penultimate) do command-> make copy" << endl;
  cout <<"##################################################################" << endl;
  cout << endl;
  return 0;
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
