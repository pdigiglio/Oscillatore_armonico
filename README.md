Descrizione
=======================

Sorgenti per lo studio numerico di un oscillatore armonico quantistico mono-dimensionale (di cui si conosce una soluzione analitica) tramite la generazione di una catena di Markov.

## Algoritmo

Il programma è un'implementazione di un [algoritmo Metropolis] [1] (di tipo Montecarlo).
L'algoritmo Metropolis è un tipo di algoritmo che soddisfa la condizione di _micro-reversibilità_ e quella di _bilancio dettagliato_.
Quella di bilancio dettagliato è condizione sufficiente per l'ergodicità dell'algoritmo.


## Scopo del programma

Il programma ha come scopo la stima di:

* Differenza tra l'energia _E_ dello stato fondamentale e quella del primo stato eccitato;
* Elemento di matrice della posizione _x_ tra lo stato energetico fondamentale _e<sub>0</sub>_ e il primo eccitato _e<sub>1</sub>_.

Medie ed incertezze su queste osservabili sono ottenute mediante il metodo Jackknife.
In in primo momento i correlatori sono suddivisi in _bin_ di dimensione tale da potersi considerare de-correlati.
Questa dimensione è dell'ordine di dieci volte il tempo di auto-correlazione dei correlatori.
Successivamente, dai bin si ricavano i _cluster_ necessari nel metodo Jackknife.


### Differenza d'energia

Per calcolare la differenza d'energia _E_, il programma calcola il __correlatore__ tra gli operatori posizione _x_ ai tempi _i_ e _j_ sul reticolo, cioè

![equation](http://bit.ly/1ilDhES)

Al denominatore compare la _funzione di partizione_ e _S_ è l'azione del sistema.
Gli integrali che compaiono nel correlatore vengono stimati tramite il metodo Montecarlo con il campionamento d'importanza.


<!--

## Features

Il programma permette di calcolare gli auto-correlatori 
---->

[1]: http://it.wikipedia.org/wiki/Algoritmo_di_Metropolis-Hastings "Algoritmo Metropolis su Wikipedia"
