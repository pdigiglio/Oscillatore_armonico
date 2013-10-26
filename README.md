Descrizione
=======================

Sorgenti per lo studio numerico di un oscillatore armonico quantistico mono-dimensionale (di cui si conosce una soluzione analitica).

## Algoritmo

Il programma è un'implementazione di un [algoritmo Metropolis] [1] (di tipo Montecarlo).
L'algoritmo Metropolis è un tipo di algoritmo che soddisfa la condizione di _micro-reversibilità_ e quella di _bilancio dettagliato_.
Quella di bilancio dettagliato è condizione sufficiente per l'ergodicità dell'algoritmo.


## Scopo del programma

Il programma ha come scopo la stima di:

* Differenza tra l'energia dello stato fondamentale e quella del primo stato eccitato;
* Elemento di matrice della posizione _x_ tra lo stato energetico fondamentale _e<sub>0</sub>_ e il primo eccitato _e<sub>1</sub>_.

Medie ed incertezze su queste osservabili sono ottenute mediante il metodo Jackknife.

<!--

## Features

Il programma permette di calcolare gli auto-correlatori 
---->


[1]: http://it.wikipedia.org/wiki/Algoritmo_di_Metropolis-Hastings "Algoritmo Metropolis su Wikipedia"

