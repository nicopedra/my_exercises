#include "funz.cpp"

using namespace std;

 int main (int argc, char *argv[]) {
//###########################################################
Random rnd = random_initialization(3);
//###########################################################

//############################# PART 1 ##################
 int M=10000; //number of total throwns
 int N=100; //number of blocks for data blocking
 vector<double> r;
 vector<double> r_mean;
 for (int i=0;i<M;i++) r.push_back(rnd.Rannyu()); //generate random numbers
 for (int i=0;i<N;i++) 
	 r_mean.push_back( mean(r, (i+1)*M/N, i*M/N )); //doing the mean in each block
 					   		//and save result
 data_blocking(N,r_mean,0.5,"datir.txt");//blocking method 

//############################# PART 2 ######################
 
 vector<double> sigma_r;
 vector<double> sigma_mean;
 for (int i=0;i<M;i++) //for each random number calculate sigma 
	 sigma_r.push_back((r[i]-0.5)*(r[i]-0.5));//0.5 is mu
 for (int i=0;i<N;i++)
	 sigma_mean.push_back( mean(sigma_r, (i+1)*M/N, i*M/N ));
 data_blocking(N,sigma_mean,1./12,"datisigma.txt"); 

//############################# PART 3  #####################

 M=100; //number of sub-intervals in [0,1]
 int it=100; //iterations number
 int n=10000; //number of throwns for each iteration
 double lint = 1./M; //lenght of each sub-interval
 vector<double> chi_tot;
 vector<int> hist(M,0);//setting all to 0
 r.clear();// I want to use again this vector

 for (int i=0;i<it;i++){ //cycle over iterations
 	for (int k=0;k<n;k++) r.push_back(rnd.Rannyu()); //generate random numbers
	for (int l=0;l<M;l++) //cycle to fill the histogram
		for (int j=0;j<n;j++) 
			if(r[j]<(l+1)*lint && r[j]>=l*lint) hist[l]++;
	chi_tot.push_back(chiquadro(hist,n/M));//save each chi-quadro in a vector
	hist.clear();
	hist.resize(M,0);//re-initialize each element of hist to 0
	r.clear();	
 }

 //saving results in file
 print_vector(chi_tot,"chiquadro.txt");
 
 system("mv dati* chiquadro.txt ../");//move the files in the jupyter-notebook directory 

 rnd.SaveSeed();//save seed in seed.out

 return 0;

}



