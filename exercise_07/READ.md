Nella cartella in cui vi trovate sono presenti tutte le cartelle e i file necessari all'analisi dati con il jupyter notebook.  
nella cartella MonteCarlo\_NVT sono presenti invece necessari al codice c++

il codice prende in input il file input.dat da cui leggere i parametri che gli servono. 
dentro tale file è specificato in che ordine mettere i parametri e cosa significano.
Da terminale bisogna dare il parametro restart:  
- posto uguale a 0 : l'eseguibile legge da config.fcc la configurazione iniziale
- posto uguale a 1 : l'eseguibile legge da config.final la configurazione della precedente simulazione  

NOTA: la funzione input() legge solo input.dat., per cui prima di iniziare una simulazione eseguire comando:  
cp input\input.cosa\_voglio\_simulare input.dat (nella cartella ci sono i file input che possono essere necessari)

Monte\_Carlo\_NVT.cpp : contiene il main
Monte\_Carlo\_NVT.h dichiarazione e implementazione di funzioni e variabili

l'eseguibile clean.sh rimuove i file generati dal codice, 
che quindi prima di essere cancellati (se utili) vanno spostati nella cartella oppurtuna per l'analisi dati

sono presenti altri file config.0 e config.final
questi sono necessari per decidere quale configurazione assegnare a r(t0).
restart = 0 ---> r(t0) = config.fcc (configurazione reticolo fcc)
restart = 1 ---> r(t0) = config.final (ultima configurazione ottenuta dalla simulazione precedente)
config.final viene sovrascritta al termine di ogni simulazione.
 
con il comando *make copy* si copia config.final in config.0, chiamato prima che config.final venga sovrascritta  
per salvare in config.0 la configurazione da cui si era partiti

OSS: prima di eseguire l'analisi si sono svolte simulazioni una di seguito all'altra per far raggiungere
al sistema la temperatura target,cioè quella a cui si vuole studiare la fisica.   
Per fare questo occorre eseguire diverse simulazione una di seguito all'altra partendo dall'ultima simulazione ottenuta  
nella simulazione precedente.

