Nella cartella ISING\_1D c'è il codice relativo all'esercizio  
della sesta esercitazione.
Nella cartella in cui vi trovate sono presenti inoltre sia il jupyter notebook
necessario per l'analisi dati e due cartelle (Gibbs e Metropolis) contenenti 
i dati necessari all'analisi 

la cartella ISING\_1D è organizzata nel seguente modo:
0)Makefile
1)Monte\_Carlo\_ISING\_1D.cpp
2)Monte\_Carlo\_ISING\_1D.h
3)simulazione.sh
4)file utili per il generatore di numeri casuali

il (0) è il makefile (basta eseguire make e compila)
(1) è il main con anche le implementazioni delle funzioni 
(2) dichiarazione delle funzioni CON LORO DESCRIZIONE. per sapere cosa fa una
certa funzione occorre cercarla in questo file con la rispettiva descrizione
(3) bash script che fa partire una simulazione per un certo range di 
temperatura

nota: bisogna dare da linea di comando la temperatura desiderata e
il parametro restart, per ripartire dalla configurazione delle simulazione
precedente 

