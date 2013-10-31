Descrizione
=======================

Sorgenti per lo studio numerico di un oscillatore armonico quantistico mono-dimensionale (di cui si conosce una soluzione analitica) tramite la generazione di una catena di Markov.
Il programma ha come scopo la stima di:

* Differenza tra l'energia _E_ dello stato fondamentale e quella del primo stato eccitato;
* Elemento di matrice della posizione _x_ tra lo stato energetico fondamentale _e<sub>0</sub>_ e il primo eccitato _e<sub>1</sub>_.


## Algoritmo

Il programma è un'implementazione di un [algoritmo Metropolis] [1] (di tipo Montecarlo).
L'algoritmo Metropolis è un tipo di algoritmo che soddisfa la condizione di _micro-reversibilità_ e quella di _bilancio dettagliato_.
Quella di bilancio dettagliato è condizione sufficiente per l'ergodicità dell'algoritmo.


Medie ed incertezze su differenza di energia _E_ ed elemento di matrice sono ottenute mediante il metodo Jackknife.
In un primo momento i correlatori sono suddivisi in _bin_ di dimensione tale da potersi considerare de-correlati.
Questa dimensione è dell'ordine di dieci volte il tempo di auto-correlazione dei correlatori.
Successivamente, dai bin si ricavano i _cluster_ necessari nel metodo Jackknife.


### Calcolo delle osservabili

Per calcolare la differenza d'energia _E_, il programma calcola il __correlatore__ tra gli operatori posizione _x_ ai tempi _i_ e _j_ sul reticolo, cioè

![correlator](http://bit.ly/1ilDhES)

Al denominatore compare la _funzione di partizione_ e _S_ è l'azione del sistema.
Gli integrali che compaiono nel correlatore vengono stimati tramite il metodo Montecarlo con il campionamento d'importanza.

<!--
La differenza di energia si ottiene dai correlatori tramite la formula:

![energy](http://www.sciweavers.org/tex2img.php?eq=%5CDelta%20E%28t%29%20%3D%20%5Cfrac%7Bc%28t%2B1%29%20%2B%20c%28t-1%29%7D%7B2c%28t%29%7D%5Cqquad%20t%3A%3D%20%7Ci-j%7C&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)


L'elemento di matrice, invece, si ottiene da

--->

[1]: http://it.wikipedia.org/wiki/Algoritmo_di_Metropolis-Hastings "Algoritmo Metropolis su Wikipedia"
