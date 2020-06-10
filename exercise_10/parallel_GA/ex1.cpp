//module load mpi/mpich-3.2-x86_64
//mpicxx ex1.cpp
//mpiexec -np 4 ./ex1.out
//plot 'picture_path0.txt' u 2:3 w l, 'picture_path3.txt' u 2:3 w l, 'picture_path2.txt' u 2:3 w l, 'picture_path1.txt' u 2:3 w l

#include "funz.cpp"

#define square
//se voglio usare solo p=1 non uso swap
#define swap

using namespace std;

 int main (int argc, char *argv[]) {
//###########################################################
//################################# EXERCISE 1 ##################
 iterations=20000;
 int cities = 32;
 Random rnd = random_initialization(2);
 int size_population = cities;
 
 N_migr = 200; 
#ifdef square
 if_square = true;
 initialize_square(rnd,cities);
#else
 if_square = false;
 initialize_circon(rnd,cities);
#endif
 
 int processes,rank;
 MPI_Init(&argc,&argv);
 MPI_Comm_size(MPI_COMM_WORLD,&processes);
 MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 MPI_Status stat;

 double tstart = MPI_Wtime();

 Population P(size_population,cities,0.2,0.7);
 int itag = 1;
 int rank_a;
 int rank_b;
 vector<int> best_a(cities);
 vector<int> best_b(cities);

 P.fill_initial_population(rank);
 P.order_chrom();

 if (rank == 0 ) {
		cout << "---------------------------------------" << endl;
#ifdef square
		cout << cities << " cities inside a square" << endl;
#else
		cout << cities << " cities placed on a circumference" << endl;
#endif
		cout << "iterations " << iterations << endl;
		cout << "N_migr " << N_migr << endl;
		cout << "population size: " << size_population << endl;
		cout << "genes per-chromosome: " << cities << endl;
		cout << "probability crossover: " << P.get_pc() << endl;
		cout << "probability mutation: " << P.get_pm() << endl;
		cout << "---------------------------------------" << endl;
	}
	
 for (int i = 0;i<iterations;i++) {

	if (i%1000 == 0 && rank == 0) {
		cout << "***************************" << endl;
		cout << "iterazione numero: " << i << endl;
		cout << "***************************" << endl;
	}
	//migration
#ifdef swap
//scambio dei migliori tra 2 nodi scelti random	
	if(i%N_migr == 0){
		//scelgo random il primo nodo (a)
		rank_a = int(rnd.Rannyu(0.,processes));
		//scelgo il secondo nodo (b)
		rank_b = int(rnd.Rannyu(0.,processes));
		//diverso da (a)
		while (rank_b == rank_a) {
			rank_b = int(rnd.Rannyu(0.,processes));
		}
		//se sono il nodo (a) salvo il mio best in best_a
		if (rank == rank_a) {
			best_a = P.get_chromo(0);
			//lo mando a (b)
			MPI_Send(&best_a[0],cities,MPI_INTEGER,rank_b,itag,MPI_COMM_WORLD);
			//mi preparo a riceve best di (b)
			MPI_Recv(&best_b[0],cities,MPI_INTEGER,rank_b,itag,MPI_COMM_WORLD,&stat);
			//ricevuto lo metto nella popolazione al posto del primo
			P.scambio(best_b);
		}
		//il nodo (b) fa la stessa cosa
		else if (rank == rank_b) {
			best_b = P.get_chromo(0);
			MPI_Send(&best_b[0],cities,MPI_INTEGER,rank_a,itag,MPI_COMM_WORLD);
			MPI_Recv(&best_a[0],cities,MPI_INTEGER,rank_a,itag,MPI_COMM_WORLD,&stat);
			P.scambio(best_a); 
		}
		P.order_chrom();
	}
#endif
	 P.print_average_path();
	 P.print_best_path();
	 //selection
	 P.selection();
	 //crossover + mutation
	 P.crossover();
	 P.update();
 }

P.order_chrom();
P.print_average_path();
P.print_best_path();
P.print_result(rank);

picture_path(P.get_chromo(0),"picture_path"+to_string(rank)+".txt");

#ifdef swap

P.best_list(P.get_chromo(0));

//ogni nodo manda al nodo 0 il suo migliore
//il nodo zero salva nei bigFour i migliori 4 e stampa il cammino del migliore
for (int i=1;i<processes;i++) {

        best_a = P.get_chromo(0);
	if (rank == i) {
                        MPI_Send(&best_a[0],cities,MPI_INTEGER,0,itag,MPI_COMM_WORLD);
                }
        else if (rank == 0) {
                       MPI_Recv(&best_a[0],cities,MPI_INTEGER,i,itag,MPI_COMM_WORLD,&stat);
		       P.best_list(best_a);
		       picture_path(P.get_theRealOne(0),"picture_final_path.txt");
	}
}
#endif

double tend = MPI_Wtime();
double dt = tend-tstart;

#ifdef square
ofstream Time("parallel_GA_square.time",ios::out);
#else
ofstream Time("parallel_GA_circumference.time",ios::out);
#endif

Time << dt << endl;
Time.close();

#ifdef swap
if (rank==0) { 
	cout << "tempo trascorso: " << dt << endl;
	cout << "best cost_function: " << cost_function(P.get_theRealOne(0) ) << endl;
}
#endif

MPI_Finalize();

rnd.SaveSeed();

 return 0;

}



