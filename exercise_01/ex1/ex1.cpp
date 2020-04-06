#include "funz.cpp"

using namespace std;

 int main (int argc, char *argv[]) {
//###########################################################
Random rnd = random_initialization();
//###########################################################

//############################# PART 1 ##################
 int M=10000; //number of total throwns
 int N=100; //number of blocks
 int L=M/N; 
 vector<double> r;
 vector<double> r_mean;
 for (int i=0;i<M;i++) r.push_back(rnd.Rannyu()); //generate random numbers
 for (int i=0;i<N;i++) 
	 r_mean.push_back( mean(r, (i+1)*L, i*L )); //doing the mean in each block
 					   		//and save result
 data_blocking(N,r_mean,0.5,"datir.txt"); 

//############################# PART 2 ######################
 
 vector<double> sigma_r;
 vector<double> sigma_mean;
 for (int i=0;i<M;i++) //for each random number calculate sigma 
	 sigma_r.push_back((r[i]-0.5)*(r[i]-0.5));
 for (int i=0;i<N;i++)
	 sigma_mean.push_back( mean(sigma_r, (i+1)*L, i*L ));
 data_blocking(N,sigma_mean,1./12,"datisigma.txt"); 

//############################# PART 3  #####################

 //N = number of sub-intervals in [0,1]
 int it=100; //iterations number
 M=10000; //number of throwns for each iteration
 double lint = 1./N; //lenght of each sub-interval
 vector<double> chi_tot;
 vector<int> hist(N,0);
 L = M/N;
 r.clear();// I want to use again this vector

 for (int i=0;i<it;i++){ //cycle over iterations
 	for (int k=0;k<M;k++) r.push_back(rnd.Rannyu()); //generate random numbers
	for (int l=0;l<N;l++) //cycle to fill the histogram
		for (int j=0;j<M;j++) 
			if(r[j]<(l+1)*lint && r[j]>=l*lint) hist[l]++;
	chi_tot.push_back(chiquadro(hist,L));//save each chi-quadro in a vector
	hist.clear();
	hist.resize(N,0);//re-initialize each element of hist to 0
	r.clear();	
 }

 print_vector(chi_tot,"chiquadro.txt");
 
 system("mv dati* chiquadro.txt ../");//move the files in the directory 
 				      //in which there's the jupyter notebook
 rnd.SaveSeed();

 return 0;

}



