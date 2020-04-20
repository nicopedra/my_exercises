nella cartella ex_arma sono presenti tutti i file necessari per il codice in cpp. 
le altre cartelle presenti contengono i file necessari al jupyter notebook per svolgere l'analisi dati. 
nella cartella video_vmd è presente una breve visualizzazione ottenuta tramite l'utilizzo del programma VMD
(visual molecular dynamics, https://www.ks.uiuc.edu/Research/vmd/) della simulazione di un liquido generico.

cartella ex_arma:
  il nome arma deriva dal fatto che ho modificato il codice originale, riscrivendolo utilizzando la libreria
  armadillo (http://arma.sourceforge.net/). Ho tenuto comunque una copia del codice scritto senza l'utilizzo di tale libreria. 
  Svolgendo l'analisi con i dati ottenuti da entrambi i codici ho verificato di non avere differenze qualitative
  tra un codice e l'altro.

il codice prende in input il file input.dat da cui leggere i parametri che gli servono. 
dentro tale file è specificato in che ordine mettere i parametri e cosa significano.
NOTA: non bisogna dare niente da terminale. sarà la funzione input (chiamata a inizio main)
a leggere il file input.dat. Per cui se si vuole modificare il parametro restart 
(leggere dentro input.dat per maggiori info) bisogna modificare input.dat
nella cartella sono presenti più file input. , ma la funzione input() legge solo input.dat., 
per cui prima di iniziare una simulazione eseguire comando: cp input.cosa_voglio_simulare input.dat

MolDyn_NVE.cpp : contiene il main
MolDyn_NVE.h dichiarazione di funzioni e variabili
funz.cpp : implementazione delle funzioni utilizzate
directory old_config: sono presenti vecchie configurazioni utili di sistemi che hanno già raggiunto
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
in questo modo una configurazione preferita viene salvata in questi file e in caso di errore non è andata
persa

OSS: prima di eseguire l'analisi si sono svolte simulazioni una di seguito all'altra per far raggiungere
al sistema la temperatura target, quella a cui si vuole studiare la fisica
per cui la simulazione 0 dovrà essere fatta utilizzando il metodo velocity random (config.fcc , restart=0)
e tutte le altri invece usando il metodo di lettura delle configurazioni della simulazione precedente (restart=1).
all'inizio del jupyter notebook vengono mostrati i risulati di queste simulazioni.
il gas è stato il più complicato
 
verrà anche generato un file traj.xyz, utile solo per visualizzare le molecole tramite vmd
si troverà commentata questa funzione nel main

al termine del main viene chiamata la funzione data_blocking. 
ci sono due possibilità di esecuzione (una delle due è commentata, sarà l'esecuture a decidere quale usare).
	-) data blocking che restituisce l'analisi dati in unità LJ (il primo, che prende solo un argomento,
	   gli altri di default vengono messi a 0)
	-) data blocking che restituisce l'analisi dati in unità Sistema Internazionale (inserendo gli altri parametri) 

altre note interessanti:
da terminale verrà stampato anche il tempo di esecuzione del codice.
confronti--> simulazione di un solido, 100000 step:
		-) codice vecchio: 350secondi
		-) codice arma: 148secondi
per provare ad utilizzare i float basta sostituire ogni double con float e ogni mat con fmat (guadagno circa 7secondi)


