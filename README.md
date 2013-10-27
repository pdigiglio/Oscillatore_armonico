Descrizione
=======================

Sorgenti per lo studio numerico di un oscillatore armonico quantistico mono-dimensionale (di cui si conosce una soluzione analitica).

## Algoritmo

Il programma è un'implementazione di un [algoritmo Metropolis] [1] (di tipo Montecarlo).
L'algoritmo Metropolis è un tipo di algoritmo che soddisfa la condizione di _micro-reversibilità_ e quella di _bilancio dettagliato_.
Quella di bilancio dettagliato è condizione sufficiente per l'ergodicità dell'algoritmo.


## Scopo del programma

Il programma ha come scopo la stima di:

* Differenza tra l'energia _E_ dello stato fondamentale e quella del primo stato eccitato;
* Elemento di matrice della posizione _x_ tra lo stato energetico fondamentale _e<sub>0</sub>_ e il primo eccitato _e<sub>1</sub>_.

Medie ed incertezze su queste osservabili sono ottenute mediante il metodo Jackknife.


### Differenza d'energia

Per calcolare la differenza d'energia _E_, il programma calcola il __correlatore__ definito come 

![equation](http://bit.ly/1ilDhES)

Il programma stima anche il correlatore _c(i,j)_


<!--
<img align="center" src="http://www.sciweavers.org/tex2img.php?eq=c%28i%2Cj%29%20%3A%3D%20%5Cdfrac%7B%20%5Cdisplaystyle%5Cint%20x_ix_j%20%5C%2C%20%5Cmathrm%7Be%7D%5E%7B-S%2F%5Chbar%7D%5Cprod_%7Bk%3D0%7D%5E%7BN-1%7D%5Cmathdm%7Bd%7D%20x_k%7D%7B%5Cdisplaystyle%20%5Cint%20%5Cmathrm%7Be%7D%5E%7B-S%2F%5Chbar%7D%5Cprod_%7Bk%3D0%7D%5E%7BN-1%7D%5Cmathdm%7Bd%7D%20x_k%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" alt="equation">
--->

<!--

## Features

Il programma permette di calcolare gli auto-correlatori 
---->

[1]: http://it.wikipedia.org/wiki/Algoritmo_di_Metropolis-Hastings "Algoritmo Metropolis su Wikipedia"
