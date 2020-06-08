FLAGS = -std=c++11 -DARMA_DONT_USE_WRAPPER -larmadillo -llapack -lblas -fopenmp -O3

compila: MolDyn_NVE.cpp funz.cpp MolDyn_NVE.h 
	g++ MolDyn_NVE.cpp -o MolDyn_NVE  $(FLAGS) 

clean:  
	rm MolDyn_NVE

copy:
	cp config.0 old.0
	cp config.final old.final
oldcopy:
	cp old.0 config.0
	cp old.final config.final	
