nella cartella ex\_arma sono presenti tutti i file necessari per il codice in cpp. 
le altre cartelle presenti contengono i file necessari al jupyter notebook per svolgere l'analisi dati. 
nella cartella video\_vmd è presente una breve visualizzazione ottenuta tramite l'utilizzo del programma VMD
(visual molecular dynamics, https://www.ks.uiuc.edu/Research/vmd/) della simulazione di un liquido generico.

cartella ex\_arma:
  il nome arma deriva dal fatto che ho modificato il codice originale, riscrivendolo utilizzando la libreria
  armadillo (http://arma.sourceforge.net/). Ho tenuto comunque una copia del codice scritto senza l'utilizzo di tale libreria. 
  Svolgendo l'analisi con i dati ottenuti da entrambi i codici ho verificato di non avere differenze qualitative e quantitative
  tra un codice e l'altro.

il codice prende in input il file input.dat da cui leggere i parametri che gli servono. 
dentro tale file è specificato in che ordine mettere i parametri e cosa significano.
Da terminale bisogna dare il parametro restart:  
- posto uguale a 0 : l'eseguibile legge da config.fcc la configurazione iniziale, e genera con velocità random le configurazioni old
- posto uguale a 1 : l'eseguibile legge da config.0 e config.final le configurazioni della precedente simulazione  

NOTA: la funzione input() legge solo input.dat., per cui prima di iniziare una simulazione eseguire comando:  
cp input\input.cosa\_voglio\_simulare input.dat (nella cartella ci sono i file input che possono essere necessari)

MolDyn\_NVE.cpp : contiene il main
MolDyn\_NVE.h dichiarazione di funzioni e variabili
funz.cpp : implementazione delle funzioni utilizzate
directory old\_config: sono presenti vecchie configurazioni utili di sistemi che hanno già raggiunto
			la temperatura target (relativa al loro file input.)

l'eseguibile clean.sh rimuove i file generati dal codice, 
che quindi prima di essere cancellati (se utili) vanno spostati nella cartella oppurtuna per l'analisi dati

sono presenti altri file config. e old.
questi sono necessari per decidere quale configurazione assegnare a r(t0+dt) (t0 = step 0 della simulazione) e r(t0).
restart = 0 ---> r(t0+dt) = config.fcc (configurazione reticolo fcc)
		 r(t0) = metodo velocità random
restart = 1 ---> r(t0+dt) = config.0 (ultima configurazione ottenuta dalla simulazione precedente)
		 r(t0) = config.final (penultima configurazione ottenuta dalla simulazione precedente)
config.0 e config.final vengono sovrascritti al termine di ogni simulazione.
 
Per cui se al termine di una simulazione ci si accorge di aver commesso un errore 
(ad esempio inserimento di parametri sbagliati) non è possibile recuperare le
configurazioni da cui si era partiti dai file config.
per questo motivo al termine della simulazione compare un messaggio di attenzione che dice:
nel caso in cui si voglia salvare le configurazioni ottenute da questa simulazione allora fare comando 
*make copy*. questo copia config.0 in old.0, e config.final in old.final

OSS: prima di eseguire l'analisi si sono svolte simulazioni una di seguito all'altra per far raggiungere
al sistema la temperatura target,cioè quella a cui si vuole studiare la fisica.   
Per fare questo all'interno del main è presente la parola: *equilibration*. Se questa è definita allora l'eseguibile  
compie un ciclo while controllando la differenza tra la temperatura target e la temperatura media raggiunta dal sistema  
(valutata con il metodo blocking). Se il modulo di questa differenza è entro l'errore calcolato con il metodo blocking allora  
il programma esce dal while, salva (con *make copy*) la configurazione da cui si è partiti per ottenere l'equilibrazione in   
old.0 e old.final. Da queste configurazione si fa partire la nuova simulazione per generare i dati di cui si vuole fare  
l'analisi  
 
verrà anche generato un file traj.xyz, utile solo per visualizzare le molecole tramite vmd
si troverà commentata questa funzione nel main

al termine del main viene chiamata la funzione data\_blocking. 

(per provare ad utilizzare i float basta sostituire ogni double con float e ogni mat con fmat (guadagno circa 7secondi))


