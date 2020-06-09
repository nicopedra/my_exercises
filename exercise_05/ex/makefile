FLAGS = -std=c++11 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas -fopenmp -O3

compila: ex1.cpp funz.cpp funz.h 
	g++ ex1.cpp -o ex1 $(FLAGS) 

clean:  
	rm ex1  
