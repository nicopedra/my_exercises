FLAGS = -std=c++17 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas -fopenmp -O3 -pg

compila: ex1.cpp funz.cpp funz.h 
	g++ ex1.cpp -o ex1 $(FLAGS) 

clean:  
	rm ex1 picture_path.txt best_path.txt path_average.txt
