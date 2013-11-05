/* C standard libraries */
#include <stdlib.h>
#include <stdio.h>

#ifndef __OSCILLATORE_H__
#define __OSCILLATORE_H__


#define MAIN_PROGRAM
#include "global.h"

class oscillatore {
	public:
		/* random initial positions */
		oscillatore ( void );
		/* all initial positions are set to 'value' */
		oscillatore ( long double value );
		/* destructor */
		virtual ~oscillatore (void);
		
		/* XXX */
		void correlator (void);
		void print_c (void);
		void Cluster (void);
		void round ( long double val, long double err );

		/* METODI */

		/* stampa lo stato per un grafico */
		void plot_state (void);
		/* stampa il correlatore per un grafico */
		void plot_correlator (void);
		/* stampa gli autocorrelatori per un grafico */
		void plot_autocorrelator (void);
		
		/* restituisce l'azione */
		long double get_action (void);
		/* restitiusce il num. di aggiornamenti ad ogni sweep */
		unsigned short int get_updated (void);
		/* restitiusce la posizione x[t] (t = tempo reticolo) */
		long double get_position (unsigned short int t);
		
		/* imposta il numero di correlatori in ogni bin */
		void set_BinSize (unsigned short int n);

		/* decide se la nuova configurazione è accettata */
		void sweep (void);
		/* 
		 * calcola l'autocorrelatore
		 * stop = num. a cui mi fermo (es. ne calcolo ~ 100)
		 * t = tempo reticolo
		 */
		void autocorrelator (unsigned short int t);
		/* calcola errore ed energia */
		void observables (bool plot);
		
	private:
		/* posizioni lungo il reticolo */
		long double x[N];
		/* variabili aggiornate ad ogni sweep */
		unsigned short int updated = 0;

		struct {
			/* mean array */
			long double mean[N] = {};
			/* error (SDOMs) array */
			long double err[N] = {};
			/* times I update correlator */
			unsigned int sweep = 0;
		} c; /* correlator */

		/* matrice dei bin */
		long double *bin[N] = {};
		/* Bin size (in unit of markovian time) */
		unsigned int bs = 100;

		long double *auc[N] = {};
		/* autocorrelator array */
		long double ac[30] = {};
		long double mean = 0;

		long double *corr[N] = {};
		unsigned short int count = 1;

		/*
		 * controllano gli estremi dei tempi per cui calcolare i
		 * cluster e le osservabili
		 */
		unsigned short int tMin = 1, tMax = 8;
		
		long double S; /* azione */

		/* potenziale (armonico) */
		long double V (long double x);
		/* differenza d'azione tra due configurazioni */
		long double diff (long double step, unsigned short int i);
};

#endif /* __OSCILLATORE_H__ */
