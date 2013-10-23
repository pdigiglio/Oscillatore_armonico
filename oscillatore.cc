#include "./oscillatore.h"

using namespace std;

	/* costruttore default */
	oscillatore::oscillatore ( void ) {
		for ( unsigned short int i = 0; i < N; i ++ ) {
			x[i] = (long double) 2 * rand() / RAND_MAX;
			c.data[i] = (long double) 0;
		}
	}
	
	/* costruttore */
	oscillatore::oscillatore ( long double value ) {
		for ( unsigned short int i = 0; i < N; i ++ ) {
			x[i] = value;
			c.data[i] = (long double) 0;
		}
	}
	
	/* distruttore */
	oscillatore::~oscillatore ( void ) {}
	

	void oscillatore::plot_state ( void ) {
		for ( unsigned short int i = 0; i < N; i ++ )
			printf("%u\t%Lf\n", i, x[i]);
	}
	
	void oscillatore::plot_correlator ( void ) {
		/* controlla se 'c.data[t]' e' riempito (se no lo riempie) */
		fill_correlator();
		/* controlla se e' normalizzato */
		normalize_correlator();
		
		for ( unsigned short int i = 0; i < N; i ++ )
			printf("%u\t%Lg\n", i, c.data[i]);
	}

	void oscillatore::plot_autocorrelator(void) {
		if ( ac.initialized ) {
			/* lo stampo normalizzato a 'ac[0]' */
			for (unsigned short int j = 0; j < ac.length; j += 1)
				printf("%u\t%Lf\n", j, ac.data[j]/ac.data[0]);
		}
	}
	

	/* calcola e restituisce il valore dell'azione */
	long double oscillatore::get_action (void) {
		S = (long double) 0;
		for (unsigned short int i = 0; i < N; i += 1)
			S += (M*powl(x[(N + i+1)%N]-x[i], (long double) 2))/(2*a) + a*V(x[i]);
			
		return S;
	}
	
	/* restituisce le variabili aggiornate ad ogni sweep */
	unsigned short int oscillatore::get_updated (void) {
		return updated;
	}

	/* restituisce 'corr.data[k][t]' */
	long double oscillatore::get_corr (unsigned int k, unsigned short int t) {
		return corr.data[k][t];
	}
	
	/* restitiusce 'ac.data[t]' */
	long double oscillatore::get_ac (unsigned short int t) {
		return ac.data[t%ac.length];
	}
	
	/* restituisce 'cb.data[k][t]' */
	long double oscillatore::get_cb (unsigned int k, unsigned short int t) {
		return cb.data[k][t];
	}
	
	long double oscillatore::get_position (unsigned short int t) {
		/* metto 't % N' per evitare errori di segmentazione */
		return x[t%N];
	}
	
	
	/* imposto il numero di correlatori in un intervallo */
	void oscillatore::set_nMax (unsigned short int n) {
		cb.nMax = n;
	}
	

	/* decide se la nuova configurazione è accettata */
	void oscillatore::sweep (void) {
		/* azzero il numero di aggiornamenti per lo sweep */
		updated = 0;
		for (unsigned short int i = 0; i < N; i += 1) {
			/* estraggo dei numeri casuali tra 0 e 1 */
			for (unsigned short int j = 0; j < 2; j += 1)
				step[j] = (long double) rand()/RAND_MAX;
				
			/* controlla se $e^{\delta S} \ge r_2$ */
			if ( expl(-diff(step[0], i)) >= step[1] ) {
				/* 'updated' controlla gli aggiornamenti in uno sweep */
				updated++;
				x[i] += (long double) DELTA*(2*step[0] - 1);
			}
		}
	}

	/* tiene traccia dei valori del correlatore c[t] ad ogni ciclo */
	void oscillatore::save_correlator (void) {	
		/*
		 * Espando la matrice 'corr.(**data)' del correlatore.
		 * NOTA: qui 'corr.rows' rappresenta (anche) il tempo Markoviano;
		 *       all'inizio 'corr.rows' vale 0.
		 */
		corr.data = (long double **) realloc(corr.data, (corr.rows + 1)*sizeof(long double *));
		if ( corr.data == NULL )
			exit(EXIT_FAILURE);
		
		/* espando ogni puntatore per N (campioni reticolo) */
		corr.data[corr.rows] = (long double *) malloc(N*sizeof(long double));
		if ( corr.data[corr.rows] == NULL )
			exit(EXIT_FAILURE);
		
		for (unsigned short int t = 0; t < N; t += 1) {
			/* azzero i correlatori al tempo markoviano 'number' */
			corr.data[corr.rows][t] = (long double) 0;
			/* media "ciclica" */
			for (unsigned short int i = 0; i < N; i += 1)
				corr.data[corr.rows][t] += x[i]*x[(i+t)%N];
			/* normalizzo */
			corr.data[corr.rows][t] = (long double) corr.data[corr.rows][t]/N;
		}
		/* aggiorno num. righe */
		corr.rows++;
		/* aggiorno 'corr.initialized' */
		if ( !corr.initialized )
			corr.initialized = true;
	}
	
	/* 
	 * calcola l'autocorrelatore
	 * stop = num. a cui mi fermo (es. ne calcolo ~ 100)
	 * t = tempo reticolo
	 */
	void oscillatore::autocorrelator (unsigned short int t, unsigned short int stop = 100) {
		/* dipende da 'corr.data[][]' */
		if ( corr.initialized ) {
			/* controlla se c[t] e' riempito (se no lo riempie) */
			fill_correlator();
		
			/* espando il puntatore 'ac' */
			ac.data = (long double *) malloc(stop*sizeof(long double));
			if ( ac.data == NULL )
				exit(EXIT_FAILURE);

			/* numero di correlatori salvati < tempo max */
			if ( c.sweep < stop ) {
				printf("Sono stati salvati pochi correlatori, procedo comunque.\n");
				stop = c.sweep;
			}

			for (unsigned short int k = 0; k < stop; k += 1) {
				ac.data[k] = (long double) 0;
				for (unsigned int i = 0; i < c.sweep - k; i += 1)
					ac.data[k] += (long double) corr.data[i][t]*corr.data[i + k][t];
				/* normalizzo */
				ac.data[k] = (long double) ac.data[k]/(c.sweep - k);
				/*sottraggo la media */
				ac.data[k] -= powl(c.data[t], (long double) 2);
			}
			/* salvo la lunghezza di 'ac' */
			ac.length = stop;
			ac.initialized = true;		
		}
		else
			printf("Bisogna salvare lo storico dei correlatori per calcolare l'autocorrelatore\n");
	}
	
	/* divide i correlatori in bin mediando su 'nMax' tempi markoviani */
	void oscillatore::correlator_bin (void) {
		/* se 'corr.data[][]' e' stato inizializzato */
		if ( corr.initialized ) {
			/* lunghezza (n. righe) della matrice */
			cb.rows = (unsigned int) corr.rows/cb.nMax;
			/* se 'corr.rows' non divisibile per 'cb.nMax' allungo di uno 'cb.data[][N]' */
			if ( corr.rows % cb.nMax != 0 )
				cb.rows++;

			/* espando 'cb.(**data)' per "righe" */
			cb.data = (long double **) malloc(cb.rows*sizeof(long double *));
			if ( cb.data == NULL )
				exit(EXIT_FAILURE);
		
			/* 'stop' e' una variabile ausiliaria */
			unsigned int stop = cb.nMax;
			for (unsigned int k = 0; k < cb.rows; k += 1) {
				/* espando 'cb.(*data[k])' per "colonne" */ 
				cb.data[k] = (long double *) malloc(N*sizeof(long double));
				if ( cb.data[k] == NULL ) 
					exit(EXIT_FAILURE);
			
				/* evito errori di segmentazione */
				if ( (k + 1)*cb.nMax > corr.rows )
					stop = corr.rows - k*cb.nMax;
	
				/* assegno ad ogni 'cb.data[k][t]' la media (su 'k') dei 'corr.data[k][t]' */
				for (unsigned short int t = 0; t < N; t += 1) {
					/* azzero le variabili */
					cb.data[k][t] = (long double) 0;
					for (unsigned int i = 0; i < stop; i += 1)
						cb.data[k][t] += (long double) corr.data[k*cb.nMax + i][t];
					/* normalizzo */
					cb.data[k][t] = cb.data[k][t]/stop;
				}
			}
			
			if ( !cb.initialized )
				cb.initialized = true;

			if ( !cb.normalized )
				cb.normalized = true;		
		}
		else
			printf("Bisogna salvare lo storico dei correlatori per calcolare i \"bin\"!\n");
	}
	
	/* crea i cluster */
	void oscillatore::create_clusters (void) {
		if ( cb.initialized ) {
			/* normalizzo 'cb.data[][]' */
			normalize_correlator_binned();
			cluster.rows = cb.rows;
		
			/* espando la matrice 'cluster.(**data)' per "righe" */
			cluster.data = (long double **) malloc(cluster.rows*sizeof(long double *));
			if ( cluster.data == NULL )
				exit(EXIT_FAILURE);
			
			for (unsigned int k = 0; k < cluster.rows; k += 1) {
				/* espando ogni riga */
				cluster.data[k] = (long double *) malloc(N*sizeof(long double));
				if ( cluster.data[k] == NULL )
					exit(EXIT_FAILURE);
			
				/* creo i cluster */
				for (unsigned int t = 0; t < N; t += 1) {
					cluster.data[k][t] = (long double) 0;
					for (unsigned int j = 0; j < cb.rows; j += 1)
						cluster.data[k][t] += cb.data[j][t];
					/* sottraggo la misurazione k-esima */
					cluster.data[k][t] -= cb.data[k][t];
					/* normalizzo */
					cluster.data[k][t] = cluster.data[k][t]/(cb.rows - 1);
				}
			}
			cluster.initialized = true;
		}
		else
			printf("Bisogna \"binnare\" i correlatori prima di creare i cluster!\n");
	}
	
	/* calcola l'errore sui correlatori */
	void oscillatore::correlator_errors (bool plot = false) {
		if ( cluster.initialized ) {
			/* controllo se il correlatore e' riempito */
			fill_correlator();
		
			/* calcolo l'errore sui correlatori */
			for (unsigned short int j = 0; j < N; j += 1) {
				c.err[j] = (long double) 0;
				for (unsigned int b = 0; b < cluster.rows; b += 1)
					c.err[j] += powl(cluster.data[b][j] - c.data[j], 2);
				/* prendo la radice quadrata */
				c.err[j] = sqrtl( (long double) c.err[j]*(cluster.rows - 1)/cluster.rows);
				/* eventualmente stampa grafico */
				if ( plot )
					printf("%u\t%Lf\t%Lf\n", j, c.data[j], c.err[j]);
			}
		}
		else
			printf("Bisogna creare i cluster per calcolare l'errore sul correlatore!\n");
	}
	
	/* calcola la differenza di energia */
	void oscillatore::energy (bool plot = false) {
		if ( cluster.initialized ) {
			e.data = (long double **) malloc(cluster.rows*sizeof(long double *));
			if ( e.data == NULL )
				exit(EXIT_FAILURE);
				
			/* calcolo il cluster $\Delta E^k(t) $ con $t \in {2, ..., 6} */
			for (unsigned int k = 0; k < cluster.rows; k += 1) {
				e.data[k] = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
				if ( e.data[k] == NULL )
					exit(EXIT_FAILURE);
				/* creo i cluster dell'energia */
				for (unsigned short int t = e.time[0]; t <= e.time[1]; t += 1)
					e.data[k][t - e.time[0]] = acoshl((cluster.data[k][(t+1)%N] + cluster.data[k][(N+t-1)%N])/(2*cluster.data[k][t%N]));
			}

			/* energia media */			
			e.value = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
			if ( e.value == NULL )
				exit(EXIT_FAILURE);

			for (unsigned short int j = 0; j <= e.time[1] - e.time[0]; j  += 1) {
				/* azzero */
				e.value[j] = (long double) 0;
				/* calcolo la media */
				for (unsigned int s = 0; s < cluster.rows; s += 1)
					e.value[j] += e.data[s][j];
				/* normalizzo */
				e.value[j] = (long double) e.value[j]/cluster.rows;
			}
			
			/* errore */	
			e.error = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
			if ( e.error == NULL )
				exit(EXIT_FAILURE);
				
			for (unsigned short int q = 0; q <= e.time[1] - e.time[0]; q += 1) {
				/* azzero */
				e.error[q] = (long double) 0;
				/* calcolo errore */
				for (unsigned int w = 0; w < cluster.rows; w += 1)
					e.error[q] += powl(e.data[w][q] - e.value[q], 2);
				/* normalizzo e prendo la radice */
				e.error[q] = sqrtl((long double) e.error[q]*(cluster.rows - 1)/cluster.rows);
				/* eventualmente stampo il grafico */
				if ( plot )
					printf("%u\t%Lg\t%Lg\n", q + e.time[0], e.value[q], e.error[q]);
			}
			e.initialized = true;
		}
		else
			printf("Bisogna creare i cluster per calcolare l'errore sul correlatore!\n");
	}
	
	/* calcolo e'elemento di matrice */
	void oscillatore::ananas (bool plot = false) {
		if ( cluster.initialized && e.initialized ) {
			exe.data = (long double **) malloc(cluster.rows*sizeof(long double *));
			if ( exe.data == NULL )
				exit(EXIT_FAILURE);
				
			/* calcolo il cluster */
			for (unsigned int k = 0; k < cluster.rows; k += 1) {
				exe.data[k] = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
				if ( exe.data[k] == NULL )
					exit(EXIT_FAILURE);
					
				/* creo i cluster */
				for (unsigned short int t = e.time[0]; t <= e.time[1]; t += 1) {
					exe.data[k][t - e.time[0]] = sqrtl(
						cluster.data[k][t]*expl(N*a*e.data[k][t - e.time[0]]/2)/
						(2*coshl(-a*(N/2 - t)*e.data[k][t - e.time[0]]))
						);
//					printf("%Lf\t%Lf\n", exe.data[k][t - e.time[0]]);
				}
			}

			/* faccio la media */
			exe.value = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
			if ( exe.value == NULL )
				exit(EXIT_FAILURE);

			for (unsigned short int j = 0; j <= e.time[1] - e.time[0]; j  += 1) {
				/* azzero */
				exe.value[j] = (long double) 0;
				/* calcolo la media */
				for (unsigned int s = 0; s < cluster.rows; s += 1)
					exe.value[j] += exe.data[s][j];
				/* normalizzo */
				exe.value[j] = (long double) exe.value[j]/cluster.rows;
			}
			
			/* errore */	
			exe.error = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
			if ( exe.error == NULL )
				exit(EXIT_FAILURE);
				
			for (unsigned short int q = 0; q <= e.time[1] - e.time[0]; q += 1) {
				/* azzero */
				exe.error[q] = (long double) 0;
				/* calcolo errore */
				for (unsigned int w = 0; w < cluster.rows; w += 1)
					exe.error[q] += powl(e.data[w][q] - e.value[q], 2);
				/* normalizzo e prendo la radice quadrata */
				exe.error[q] = sqrtl((long double) e.error[q]*(cluster.rows - 1)/cluster.rows);
				if ( plot )
					printf("%u\t%Lg\t%Lg\n", q + e.time[0], exe.value[q], exe.error[q]);
			}
		}
	}
	
	
	/* riempie il vettore c[N] */
	void oscillatore::update_correlator (void) {		
		for (unsigned short int t = 0; t < N; t += 1) {
			/* media "ciclica" */
			for (unsigned short int i = 0; i < N; i += 1)
				c.data[t] += x[i]*x[(i+t)%N]/N;
		}
		/* il correlatore e' riempito */
		if ( !c.filled )
			c.filled = true;
		/* aggiorno il num. di volte che registro c[N] */
		c.sweep++;
	}
	
	/* crea i "bin" man mano che faccio evolvere il sistema */
	void oscillatore::update_correlator_binned (void) {
		/* 'cb.sweep' e' inizializzato a zero */
		cb.rows = (unsigned int) cb.sweep/cb.nMax + 1;
		unsigned int k = cb.rows - 1;
		
		/* ogni 'cb.nMax' sweep allungo la matrice */
		if ( cb.sweep % cb.nMax == 0 ) {
			/* espando 'cb.(**data)' per "righe" */
			cb.data = (long double **) realloc(cb.data, cb.rows*sizeof(long double *));
			if ( cb.data == NULL )
				exit(EXIT_FAILURE);

			/* espando 'cb.(**data)' per "colonne" */
			cb.data[k] = (long double *) malloc(N*sizeof(long double));
			if ( cb.data[k] == NULL ) 
				exit(EXIT_FAILURE);

			/* azzero la matrice */
			for (unsigned short int j = 0; j < N; j += 1)
				cb.data[k][j] = (long double) 0;
		}
		

		/* aggiorno 'cb.data[k][]' con i correlatori */
		for (unsigned short int t = 0; t < N; t += 1) {
			 /* media "ciclica" */
			for (unsigned short int i = 0; i < N; i += 1)
				cb.data[k][t] += x[i]*x[(i+t)%N]/N;
		}

		if ( !cb.initialized )
			cb.initialized = true;
		/* incremento numero sweep */
		cb.sweep++;
	}
	
	/* normalizza il correlatore (se non lo e') */
	void oscillatore::normalize_correlator (void) {
		if ( !c.normalized && c.filled ) {
			for (unsigned short int t = 0; t < N; t += 1)
				c.data[t] = c.data[t]/c.sweep;
			c.normalized = true;
		}
	}
	
	/* normalizza i bin */
	void oscillatore::normalize_correlator_binned (void) {
		if ( !cb.normalized && cb.initialized) {
			unsigned short int norm = cb.nMax;
			for (unsigned int k = 0; k < cb.rows; k += 1) {
				/* calcolo il coefficiente di normalizzazione per l'ultima riga */
				if ( (k + 1)*cb.nMax > cb.sweep )
					norm = cb.sweep - k*cb.nMax;
				/* normalizzo */
				for (unsigned short int t = 0; t < N; t += 1)
					cb.data[k][t] = (long double) cb.data[k][t]/norm;
			}
			cb.normalized = true;
		}
	}


	/* potenziale (armonico) */
	long double oscillatore::V (long double x) {
		/* $m\omega^2 x^2/2$ */
		return (long double) M*powl(W*x, (double) 2)/2;
//		return 0; /* particella libera */
	}
	
	/*
	 * calcola la differenza di azione tra due configurazioni 
	 * che differiscono soltanto per il termine i-esimo
	 */
	long double oscillatore::diff (long double step, unsigned short int i) {
		/*
		 * NOTA: le espressioni (N + i - 1)%N e (i+1)%N
		 * servono per assicurare che gli indici del vettore
		 * siano tutti minori di N e ciclici
		 */
		long double delta = DELTA*(2*step - 1);
		return M*delta*(2*x[i] - x[(N + i-1)%N]- x[(i+1)%N] + delta)/a + a*(V(x[i]+delta)-V(x[i]));
		/*
		 * Lascio il potenziale indicato (non lo svolgo esplicitamente
		 * come ho fatto per l'azione) in modo da poterlo cambiare per 
		 * rescrivere altri sistemi con questo algoritmo
		 */
	}
	
	/* riempie il correlatore */
	void oscillatore::fill_correlator (void) {
		/* controllo se  e' stato riempito */
		if ( !c.filled ) {
			/* se i 'cb.data[][]' esiste, uso quello */
			if ( cb.initialized ) {
				/* (eventualmente) normalizzo i bin */
				normalize_correlator_binned();
				/* calcolo i correlatori */
				for (unsigned short int t = 0; t < N; t += 1) {
					for (unsigned int k = 0; k < cb.rows; k += 1)
						c.data[t] += cb.data[k][t];
					/* normalizzo */
					c.data[t] = (long double) c.data[t]/cb.rows;
				}
			}
			/* altrimenti uso lo storico dei correlatori 'corr.data[][]' */
			else if ( corr.initialized ) {
				for (unsigned short int t = 0; t < N; t += 1) {
					for (unsigned int k = 0; k < corr.rows; k += 1)			
						c.data[t] += corr.data[k][t];
					/* normalizzo */
					c.data[t] = c.data[t]/corr.rows;
				}
				/* assegno a 'c.sweep' il numero di sweep */
				c.sweep = corr.rows;
			}
			/* se lo riempio in questa funzione e' anche normalizzato */
			c.normalized = true;
		}
	}
