/* C standard libraries */
#include <stdlib.h>
#include <stdio.h>

#ifndef __OSCILLATORE_H__
#define __OSCILLATORE_H__


#define MAIN_PROGRAM
#include "global.h"

class oscillatore {
	public:

		long double cacca[30] = {};
		/* inizializza la posizione in modo casuale */
		oscillatore (void);
		/* inizializza la posizione assegnando 'value' ad ogni componente */
		oscillatore (long double value);
		/* distruttore */
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
//		void plot_autocorrelator (void);
		
		/* restituisce l'azione */
		long double get_action (void);
		/* restitiusce il num. di aggiornamenti ad ogni sweep */
		unsigned short int get_updated (void);
		/* restituisce 'corr.data[k][t]' */
//		long double get_corr (unsigned int k, unsigned short int t);
		/* restitiusce 'ac.data[t]' */
//		long double get_ac (unsigned short int t);
		/* restituisce 'cb.data[k][t]' */
//		long double get_cb (unsigned int k, unsigned short int t);
		/* restitiusce la posizione x[t] (t = tempo reticolo) */
		long double get_position (unsigned short int t);
		
		/* imposta il numero di correlatori in ogni bin */
		void set_BinSize (unsigned short int n);

		/* decide se la nuova configurazione è accettata */
		void sweep (void);
		/* tiene traccia del correlatore c[t] ad ogni stato */
//		void save_correlator (void);
		/* 
		 * calcola l'autocorrelatore
		 * stop = num. a cui mi fermo (es. ne calcolo ~ 100)
		 * t = tempo reticolo
		 */
		void autocorrelator (unsigned short int t);
		/* divide i correlatori in bin da 'nMax' */
//		void correlator_bin (void);
		/* crea una matrice di cluster per il metodo Jacknife */
//		void create_clusters (void);
		/* calcolo gli errori sul correlatore */
//		void correlator_errors (bool plot);
		/* calcola errore ed energia */
		void observables (bool plot);
		/* calcola elemento matrice e errori */
		void matrix (bool plot);
		
		/* NOTA:
		 *  le funzioni di tipo 'update_*' richiedono che venga dato 
		 *  'normalize_*' alla fine del ciclo di sweep.
		 */
		 
		/* aggiorna il correlatore */
//		void update_correlator (void);
		/* crea i "bin" man mano che faccio evolvere il sistema */
//		void update_correlator_binned (void);
		/* normalizza il corelatore */
//		void normalize_correlator (void);
		/* normalizza i bin */
//		void normalize_correlator_binned (void);

	private:
		/* posizioni lungo il reticolo */
		long double x[N];
		/* variabili aggiornate ad ogni sweep */
		unsigned short int updated = 0;

		struct {
			/* vettore dati, vettore errori */
			long double mean[N] = {}, err[N] = {};
			/* matrice dei bin */
			long double *bin[N] = {};
			unsigned int nB = 100;
			/* dice se e'riempito o no */
//			bool filled = false;
//			/* controlla se e' normalizzato o no */
//			bool normalized = false;
			/* numero di volte che lo aggiorno */
			unsigned int sweep = 0;
		} c; /* correlatore (va diviso per 'c.sweep') */

		long double *auc[N] = {};


		struct {
			/* matriciozza */
			long double **data = NULL;
			/* inizializzata */
			bool initialized = false;
			/* lunghezza (tempo markoviano) della matriciozza */
			unsigned int rows = 0;
		} corr; /* storico dei correlatori */

		struct {
			/* vettore dati */
			long double *data;
			/* inizializzato */
			bool initialized = false;
			/* lunghezza di 'ac[]' */
			unsigned short int length;
		} ac; /* autocorrelatore */
		
		struct {
			/* matriciozza */
			long double **data = NULL;
			/* inizializzata, normalizzata */
			bool initialized = false, normalized = false;
			/* num. correlatori in ogni bin */
			unsigned short int nMax = 100;
			/* lunghezza (tempo markoviano) della matriciozza */
			unsigned int sweep = 0, rows = 0;
		} cb; /* correlatori "binnati" */ 
		
		struct {
			/* matriciozza */
			long double **data;
			/* righe della matriciozza */
			unsigned int rows;
			/* inizializzata */
			bool initialized = false;
		} cluster; /* per il metodo Jacknife */
		
		struct {
			/* cluster per l'energia */
			long double **data;
			/* limiti in cui calcolo l'energia ($t\in{2, ..., 6}$) */
			unsigned short int time[2] = {1, 15};
			/* valore e errore dell'energia ($\Delta E$)*/
			long double *value, *error;
			/* controlla se e' inizializzata */
			bool initialized = false;
		} e; /* differenza energia */

		struct {
			/* cluster per elemento matrice */
			long double **data;
			/* valore e errore elemento matrice */
			long double *value, *error;
		} exe; /* elemento di matrice $|<E_0|x|e_1>|$ */
		
		long double S; /* azione */

		/* potenziale (armonico) */
		long double V (long double x);
		/* differenza d'azione tra due configurazioni */
		long double diff (long double step, unsigned short int i);
		/* riempie il correlatore (se non e' gia' pieno) */
		void fill_correlator (void);
};

#endif /* __OSCILLATORE_H__ */
