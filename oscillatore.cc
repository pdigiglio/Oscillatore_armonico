/* C standard math library */
#include <math.h>

/* Class header file */
#include "oscillatore.h"

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: oscillatore
 * Description: [default ctor] positions initialized at random values.
 * 				Here memory for the first correlators is allocated.
 * ------------------------------------------------------------------
 */
oscillatore::oscillatore ( void ) {
//	#pragma omp parallel for
	for ( unsigned short int t = 0; t < N; t ++ ) {
		*( x + t ) = (long double) 2 * rand() / RAND_MAX;

		/* alloco la memoria per i primi N bin del correlatore */
		*( c.bin + t ) = (long double *) malloc( sizeof(long double) );
		if ( *( c.bin + t ) == NULL ) {
			fprintf ( stderr, "[default ctor] dynamic memory allocation failed\n" );
			exit(EXIT_FAILURE);
		}
	}
} /* -----  end of method oscillatore::oscillatore (def. ctor)  ----- */

	
/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: oscillatore
 * Description: [ctor] all positions initialized at value 'value'
 * ------------------------------------------------------------------
 */
oscillatore::oscillatore ( long double value ) {
//	#pragma omp parallel for
	for ( unsigned short int t = 0; t < N; t ++ ) {
		*( x + t ) = value;

		/* alloco la memoria per i primi N bin del correlatore */
		*( c.bin + t ) = (long double *) malloc( sizeof(long double) );
		if ( *( c.bin + t ) == NULL ) {
			fprintf ( stderr, "[default ctor] dynamic memory allocation failed\n" );
			exit(EXIT_FAILURE);
		}
	}
} /* -----  end of method oscillatore::oscillatore (ctor)  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: ~oscillatore
 * Description: destructor
 * ------------------------------------------------------------------
 */
oscillatore::~oscillatore ( void ) {
} /* -----  end of method oscillatore::~oscillatore (dtor)  ----- */


/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: plot_state
 * Description: stampa a schermo lo stato (tutte le posizioni ai tem-
 * 				pi fisici del reticolo) del sistema
 * ------------------------------------------------------------------
 */
void
oscillatore::plot_state ( void ) {
	for ( unsigned short int i = 0; i < N; i ++ )
		printf( "%u\t%Lf\n", i, *(x + i) );
} /* -----  end of method oscillatore::plot_state  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: plot_correlator
 * Description: 
 * ------------------------------------------------------------------
 */
//void
//oscillatore::plot_correlator ( void ) {
//	/* controlla se 'c.mean[t]' e' riempito (se no lo riempie) */
//	fill_correlator();
//	/* controlla se e' normalizzato */
//	normalize_correlator();
//	
//	for ( unsigned short int i = 0; i < N; i ++ )
//		printf("%u\t%Lg\n", i, c.mean[i]);
//} /* -----  end of method oscillatore::plot_correlator  ----- */

//	void oscillatore::plot_autocorrelator(void) {
//		if ( ac.initialized ) {
//			/* lo stampo normalizzato a 'ac[0]' */
//			for (unsigned short int j = 0; j < ac.length; j += 1)
//				printf("%u\t%Lf\n", j, ac.data[j]/ac.data[0]);
//		}
//	}

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: get_action
 * Description: computes and returns action value
 * ------------------------------------------------------------------
 */
long double
oscillatore::get_action ( void ) {
	S = (long double) 0;

	/* temporary variable */
	long double tmp;
	for ( unsigned short int i = 0; i < N; i ++ ) {
		tmp = *( x + ( N + i + 1 ) % N ) - *( x + i );
		/* parte cinetica */
		S += .5 * M * tmp * tmp / a;
		/* aggiungo il potenziale */
		S += a * V( *( x + i ) );
	}

	return S;
} /* -----  end of method oscillatore::get_action  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: get_updated
 * Description: restituisce le variabili aggiornate ad ogni sweep 
 * ------------------------------------------------------------------
 */
unsigned short int
oscillatore::get_updated ( void ) {
	return updated;
} /* -----  end of method oscillatore::get_updated  ----- */

//	/* restituisce 'corr.data[k][t]' */
//	long double oscillatore::get_corr (unsigned int k, unsigned short int t) {
//		return corr.data[k][t];
//	}
//	
//	/* restitiusce 'ac.data[t]' */
//	long double oscillatore::get_ac (unsigned short int t) {
//		return ac.data[t%ac.length];
//	}
//	
//	/* restituisce 'cb.data[k][t]' */
//	long double oscillatore::get_cb (unsigned int k, unsigned short int t) {
//		return cb.data[k][t];
//	}


/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: get_position
 * Description: returns t-th position x[t]
 * ------------------------------------------------------------------
 */
long double
oscillatore::get_position ( unsigned short int t ) {
	/* metto 't % N' per evitare errori di segmentazione */
	return *( x + ( t % N ) );
} /* -----  end of method oscillatore::get_position  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: set_BinSize
 * Description: sets number of measures for each bin
 * ------------------------------------------------------------------
 */
void
oscillatore::set_BinSize ( unsigned short int size ) {
	c.nB = size;	
} /* -----  end of method oscillatore::set_BinSize  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: sweep
 * Description: updates system positions
 * ------------------------------------------------------------------
 */
void
oscillatore::sweep ( void ) {
	/* numero casuale */
	long double r;
	/* azzero il numero di aggiornamenti per lo sweep */
//	updated = 0;

	/* calcolo gli aggiornamenti */
	for ( unsigned short int i = 0; i < N; i ++ ) {
		/* estraggo un numero casuale tra 0 e 1 */
		r = (long double) rand() / RAND_MAX;
			
		/* controlla se $e^{\delta S} \ge r_2$ */
		if ( expl( - diff( r, i ) ) >= (long double) rand() / RAND_MAX ) {
			/* 'updated' controlla gli aggiornamenti in uno sweep */
//			updated ++;
			
			/* aggiorno la posizione (i + 1)-esima */
			*( x + i ) += (long double) DELTA * ( 2 * r - 1);
		}
	}
} /* -----  end of method oscillatore::sweep  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: correlator
 * Description: computes correlators while dividing them in bins. Pay
 * 				attention: all bins remain non-normalized.
 * ------------------------------------------------------------------
 */
void
oscillatore::correlator ( void ) {
	/* variabile ausiliaria che indica il numero di bin */
	unsigned int b = (unsigned) c.sweep / c.nB;

	if ( c.sweep < 30 )
		for ( unsigned short int t = 0; t < N; t ++ ) {
			/* alloco la memoria */
			*( auc + t ) = (long double *) realloc ( *(auc + t ), (c.sweep + 1) * sizeof( long double ) );
			if ( *( auc + t ) == NULL ) {
				fprintf ( stderr, "\ndynamic memory allocation failed\n" );
				exit (EXIT_FAILURE);
			}
		}

	/* calcolo i correlatori */
	for ( unsigned short int t = 0; t < N; t ++ ) {
		/* azzero la variabile ausiliaria */
		*( *( auc + t) + c.sweep % 30 ) = (long double) 0;
		/* calcolo il t-esimo correlatore */
		for ( unsigned short int i = 0; i < N; i ++ )
			*( *( auc + t) + c.sweep % 30 ) += *( x + i ) * *( x + ( i + t ) % N );

		/* aggiorno il bin b-esimo del correlatore t-esimo */
		*( *( auc + t ) + c.sweep % 30 ) = *( *( auc + t ) + c.sweep % 30 ) / N;
		*( *( c.bin + t ) + b ) += *( *( auc + t ) + c.sweep % 30 );
	}
	
	/*
	 * Aggiorno il numero di sweep e controllo se devo allungare la
	 * matrice dei bin. Questa matrice viene allungata quando tutti
	 * gli elementi sono stati sistemati in un bin
	 */
	if ( !( ++ c.sweep % c.nB ) ) {
		/* ricalcolo il valore di 'b' */
		b = (unsigned) c.sweep / c.nB;
		/* alloco la memoria per i nuovi bin */
		for ( unsigned short int t = 0; t < N; t ++ ) {
			*( c.bin + t ) = (long double *) realloc( *( c.bin + t), ( b + 1 ) * sizeof(long double) );
			/* controllo che non ci siano problemi di memoria */
			if ( c.bin + t == NULL ) {
				fprintf ( stderr, "[correlator] Dynamic memory reallocation failed\n" );
				exit(EXIT_FAILURE);
			}
			
			/* azzero le variabili (bin successivi) */
			*( *( c.bin + t ) + b ) = (long double) 0;

			/* aggiorno la media t-esima */
			*( c.mean + t ) += *( *( c.bin + t ) + b - 1 );
			/* aggiorno l'errore */
			*( c.err + t ) += powl( *( *( c.bin + t ) + b - 1 ), 2. );
		}
	}
} /* -----  end of method oscillatore::correlator  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: print_c
 * Description: 
 * ------------------------------------------------------------------
 */
void
oscillatore::print_c ( void ) {
	/* numero di bin */
	unsigned n = c.sweep / c.nB;
	fprintf( stderr, "%u %u: %u\n", c.sweep, c.nB, n );
	fprintf( stderr, "%Lg\t%Lg\n", c.bin[0][n - 1], c.bin[0][n - 2] );

	/* normalizzo medie ed errori */
	for ( unsigned short int t = 0; t < N; t ++ ) {
		/* normalizzo la media */
		*( c.mean + t ) = *( c.mean + t ) / n;

		/* calcolo l'errore */
		*( c.err + t ) = *( c.err + t ) / n ;
		*( c.err + t ) -= powl( *( c.mean + t ), 2. );
		/* uso la correzione di Bessel per l'errore */
		*( c.err + t ) = sqrtl( *( c.err + t ) / ( n - 1 ) );

		printf( "%hu\t%Lg\t%Lg\n", t, *( c.mean + t ), *( c.err + t ) );
	}

	fprintf( stderr, "N. bin: %u\n", n );
} /* -----  end of method oscillatore::print_c  ----- */


/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: cluster
 * Description: creates clusters from correlators bin
 * ------------------------------------------------------------------
 */
void
oscillatore::Cluster ( void ) {
	/* numero di bin */
	unsigned int bns = (unsigned) c.sweep / c.nB;

	/* assegno alle variabili contenenti i bin i cluster rispettivi */
	for ( unsigned int n = 0; n < bns; n ++ )
		for ( unsigned short int t = 0; t < N; t ++ ) {
			*( *( c.bin + t ) + n ) -= *( c.mean + t );
			*( *( c.bin + t ) + n ) = - *( *( c.bin + t ) + n ) / ( bns - 1 );
		}
} /* -----  end of method oscillatore::cluster  ----- */

	
	/* 
	 * calcola l'autocorrelatore
	 * stop = num. a cui mi fermo (es. ne calcolo ~ 100)
	 * t = tempo reticolo
	 */
	void oscillatore::autocorrelator ( unsigned short int t ) {
		cacca[0] = 0;
		for ( unsigned short int i = 0; i < MIN( c.sweep, 30 ); i ++ ) {
			cacca[i] += auc[0][ c.sweep % 30 ] * auc[0][ ( c.sweep + 30 - i ) % 30 ] ;
		}

//		/* dipende da 'corr.data[][]' */
//		if ( corr.initialized ) {
//			/* controlla se c[t] e' riempito (se no lo riempie) */
//			fill_correlator();
//		
//			/* espando il puntatore 'ac' */
//			ac.data = (long double *) malloc(stop*sizeof(long double));
//			if ( ac.data == NULL )
//				exit(EXIT_FAILURE);
//
//			/* numero di correlatori salvati < tempo max */
//			if ( c.sweep < stop ) {
//				printf("Sono stati salvati pochi correlatori, procedo comunque.\n");
//				stop = c.sweep;
//			}
//
//			for (unsigned short int k = 0; k < stop; k += 1) {
//				ac.data[k] = (long double) 0;
//				for (unsigned int i = 0; i < c.sweep - k; i += 1)
//					ac.data[k] += (long double) corr.data[i][t]*corr.data[i + k][t];
//				/* normalizzo */
//				ac.data[k] = (long double) ac.data[k]/(c.sweep - k);
//				/*sottraggo la media */
//				ac.data[k] -= powl(c.mean[t], (long double) 2);
//			}
//			/* salvo la lunghezza di 'ac' */
//			ac.length = stop;
//			ac.initialized = true;		
//		}
//		else
//			printf("Bisogna salvare lo storico dei correlatori per calcolare l'autocorrelatore\n");
	}
	
	/* divide i correlatori in bin mediando su 'nMax' tempi markoviani */
//	void oscillatore::correlator_bin (void) {
//		/* se 'corr.data[][]' e' stato inizializzato */
//		if ( corr.initialized ) {
//			/* lunghezza (n. righe) della matrice */
//			cb.rows = (unsigned int) corr.rows/cb.nMax;
//			/* se 'corr.rows' non divisibile per 'cb.nMax' allungo di uno 'cb.data[][N]' */
//			if ( corr.rows % cb.nMax != 0 )
//				cb.rows++;
//
//			/* espando 'cb.(**data)' per "righe" */
//			cb.data = (long double **) malloc(cb.rows*sizeof(long double *));
//			if ( cb.data == NULL )
//				exit(EXIT_FAILURE);
//		
//			/* 'stop' e' una variabile ausiliaria */
//			unsigned int stop = cb.nMax;
//			for (unsigned int k = 0; k < cb.rows; k += 1) {
//				/* espando 'cb.(*data[k])' per "colonne" */ 
//				cb.data[k] = (long double *) malloc(N*sizeof(long double));
//				if ( cb.data[k] == NULL ) 
//					exit(EXIT_FAILURE);
//			
//				/* evito errori di segmentazione */
//				if ( (k + 1)*cb.nMax > corr.rows )
//					stop = corr.rows - k*cb.nMax;
//	
//				/* assegno ad ogni 'cb.data[k][t]' la media (su 'k') dei 'corr.data[k][t]' */
//				for (unsigned short int t = 0; t < N; t += 1) {
//					/* azzero le variabili */
//					cb.data[k][t] = (long double) 0;
//					for (unsigned int i = 0; i < stop; i += 1)
//						cb.data[k][t] += (long double) corr.data[k*cb.nMax + i][t];
//					/* normalizzo */
//					cb.data[k][t] = cb.data[k][t]/stop;
//				}
//			}
//			
//			if ( !cb.initialized )
//				cb.initialized = true;
//
//			if ( !cb.normalized )
//				cb.normalized = true;		
//		}
//		else
//			printf("Bisogna salvare lo storico dei correlatori per calcolare i \"bin\"!\n");
//	}
	
	/* crea i cluster */
//	void oscillatore::create_clusters (void) {
//		if ( cb.initialized ) {
//			/* normalizzo 'cb.data[][]' */
//			normalize_correlator_binned();
//			cluster.rows = cb.rows;
//		
//			/* espando la matrice 'cluster.(**data)' per "righe" */
//			cluster.data = (long double **) malloc(cluster.rows*sizeof(long double *));
//			if ( cluster.data == NULL )
//				exit(EXIT_FAILURE);
//			
//			for (unsigned int k = 0; k < cluster.rows; k += 1) {
//				/* espando ogni riga */
//				cluster.data[k] = (long double *) malloc(N*sizeof(long double));
//				if ( cluster.data[k] == NULL )
//					exit(EXIT_FAILURE);
//			
//				/* creo i cluster */
//				for (unsigned int t = 0; t < N; t += 1) {
//					cluster.data[k][t] = (long double) 0;
//					for (unsigned int j = 0; j < cb.rows; j += 1)
//						cluster.data[k][t] += cb.data[j][t];
//					/* sottraggo la misurazione k-esima */
//					cluster.data[k][t] -= cb.data[k][t];
//					/* normalizzo */
//					cluster.data[k][t] = cluster.data[k][t]/(cb.rows - 1);
//				}
//			}
//			cluster.initialized = true;
//		}
//		else
//			printf("Bisogna \"binnare\" i correlatori prima di creare i cluster!\n");
//	}
	
	/* calcola l'errore sui correlatori */
//	void oscillatore::correlator_errors (bool plot = false) {
//		if ( cluster.initialized ) {
//			/* controllo se il correlatore e' riempito */
//			fill_correlator();
//		
//			/* calcolo l'errore sui correlatori */
//			for (unsigned short int j = 0; j < N; j += 1) {
//				c.err[j] = (long double) 0;
//				for (unsigned int b = 0; b < cluster.rows; b += 1)
//					c.err[j] += powl(cluster.data[b][j] - c.mean[j], 2);
//				/* prendo la radice quadrata */
//				c.err[j] = sqrtl( (long double) c.err[j]*(cluster.rows - 1)/cluster.rows);
//				/* eventualmente stampa grafico */
//				if ( plot )
//					printf("%u\t%Lf\t%Lf\n", j, c.mean[j], c.err[j]);
//			}
//		}
//		else
//			printf("Bisogna creare i cluster per calcolare l'errore sul correlatore!\n");
//	}
	
/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: energy
 * Description: uses clusters to obtain energy density
 * ------------------------------------------------------------------
 */
void
oscillatore::energy ( bool plot = true ) {
	/* variabile temporanea ausiliare */
	long double tmp;
	/* altre due variabili temporanee */
	long double mean, err;

	/* numero di bin */
	unsigned int bns = (unsigned) c.sweep / c.nB;

	/* calcolo le energie */
	for ( unsigned short int t = 1; t < 8; t ++ ) {
		/* azzero le varabili temporanee */
		mean = (long double) 0;
		err = (long double) 0;

		for ( unsigned int n = 0; n < bns; n ++ ) {
			/* assegno la variabile temporanea */
			tmp = *( *( c.bin + t + 1 ) + n ) + *( *( c.bin + t - 1 ) + n );
			tmp = acoshl( .5 * tmp / *( *( c.bin + t ) + n ) );

			/* aggiorno media ed errore */
			mean += tmp;
			err += powl( tmp, 2. );
		}

		/* normalizzo la media */
		mean = mean / bns;

		/* calcolo gli errori */
		err = err / bns - powl( mean, 2. );
		err = sqrtl( err * (bns - 1) );

		/* stampo a schermo i valori ottenuti in funzione di 't' */
		if( plot ) {
			printf( "%hu\t", t );
			oscillatore::round( mean, err );
			printf( "\n" );
		}
	}
} /* -----  end of method oscillatore::energy  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: matrix
 * Description: 
 * ------------------------------------------------------------------
 */
void
oscillatore::matrix ( bool plot = true ) {
	/* variabili temporanee ausiliarie */
	long double tmp, ene;
	/* altre due variabili temporanee */
	long double mean, err;

	/* numero di bin */
	unsigned int bns = (unsigned) c.sweep / c.nB;

	/* calcolo le energie */
	for ( unsigned short int t = 1; t < 8; t ++ ) {
		/* azzero le varabili temporanee */
		mean = (long double) 0;
		err = (long double) 0;

		for ( unsigned int n = 0; n < bns; n ++ ) {
			/* calcolo l'energia */
			ene = *( *( c.bin + t + 1 ) + n ) + *( *( c.bin + t - 1 ) + n );
			ene = acoshl( .5 * ene / *( *( c.bin + t ) + n ) );
			/* assegno la variabile temporanea */
			tmp = expl( - (long double) t * ene );
			tmp = *( *( c.bin + t ) + n ) / ( tmp + expl( - a * N * ene ) / tmp );

			/* aggiorno media ed errore */
			mean += sqrtl( tmp );
			err += tmp;
		}

		/* normalizzo la media */
		mean = mean / bns;

		/* calcolo gli errori */
		err = err / bns - powl( mean, 2. );
		err = sqrtl( err * ( bns - 1) );

		/* stampo a schermo i valori ottenuti in funzione di 't' */
		if( plot ) {
			printf( "%hu\t", t );
			oscillatore::round( mean, err );
			printf( "\n" );
		}
	}
} /* -----  end of method oscillatore::matrix  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: round
 * Description: 
 * ------------------------------------------------------------------
 */
void
oscillatore::round ( long double val, long double err ) {
	if (  isnan( err ) ) {
		printf( "%Lf\t%Lf", val, err );
		return;
	}

	/* dichiaro una variabile per l'esponente */
	short int exp = (short) log10( fabs( err ) );

	/* controllo che l'approssimazione sia corretta */
	while ( !( err / pow( 10, exp ) >= 1 && err / pow( 10, exp ) < 10 ) ) {
		if ( err / pow( 10, exp ) <= 1 ) exp --;
		else exp ++;
	}

	/* controllo le cifre decimali da tenere */
	if ( err / pow( 10, exp ) < 3 )
		exp --;

	printf( "%Lf\t%Lf",
				floorl( val / pow(10., exp) + 0.5) * pow(10., exp),
				floorl( err / pow(10., exp) + 0.5) * pow(10., exp)
			);
} /* -----  end of method oscillatore::round  ----- */

//	/* calcolo e'elemento di matrice */
//	void oscillatore::ananas (bool plot = false) {
//		if ( cluster.initialized && e.initialized ) {
//			exe.data = (long double **) malloc(cluster.rows*sizeof(long double *));
//			if ( exe.data == NULL )
//				exit(EXIT_FAILURE);
//				
//			/* calcolo il cluster */
//			for (unsigned int k = 0; k < cluster.rows; k += 1) {
//				exe.data[k] = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
//				if ( exe.data[k] == NULL )
//					exit(EXIT_FAILURE);
//					
//				/* creo i cluster */
//				for (unsigned short int t = e.time[0]; t <= e.time[1]; t += 1) {
//					exe.data[k][t - e.time[0]] = sqrtl(
//						cluster.data[k][t]*expl(N*a*e.data[k][t - e.time[0]]/2)/
//						(2*coshl(-a*(N/2 - t)*e.data[k][t - e.time[0]]))
//						);
////					printf("%Lf\t%Lf\n", exe.data[k][t - e.time[0]]);
//				}
//			}
//
//			/* faccio la media */
//			exe.value = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
//			if ( exe.value == NULL )
//				exit(EXIT_FAILURE);
//
//			for (unsigned short int j = 0; j <= e.time[1] - e.time[0]; j  += 1) {
//				/* azzero */
//				exe.value[j] = (long double) 0;
//				/* calcolo la media */
//				for (unsigned int s = 0; s < cluster.rows; s += 1)
//					exe.value[j] += exe.data[s][j];
//				/* normalizzo */
//				exe.value[j] = (long double) exe.value[j]/cluster.rows;
//			}
//			
//			/* errore */	
//			exe.error = (long double *) malloc((e.time[1] - e.time[0] + 1)*sizeof(long double));
//			if ( exe.error == NULL )
//				exit(EXIT_FAILURE);
//				
//			for (unsigned short int q = 0; q <= e.time[1] - e.time[0]; q += 1) {
//				/* azzero */
//				exe.error[q] = (long double) 0;
//				/* calcolo errore */
//				for (unsigned int w = 0; w < cluster.rows; w += 1)
//					exe.error[q] += powl(e.data[w][q] - e.value[q], 2);
//				/* normalizzo e prendo la radice quadrata */
//				exe.error[q] = sqrtl((long double) e.error[q]*(cluster.rows - 1)/cluster.rows);
//				if ( plot )
//					printf("%u\t%Lg\t%Lg\n", q + e.time[0], exe.value[q], exe.error[q]);
//			}
//		}
//	}
	
	
	/* riempie il vettore c[N] */
//	void oscillatore::update_correlator (void) {
//
//		for ( unsigned short int t = 0; t < N; t ++ ) {
//			/* media "ciclica" */
//			for ( unsigned short int i = 0; i < N; i ++ )
//				c.mean[t] += *( x + i ) * *( x + ( i + t ) % N );
//			
//			/* normalizzo la media */
//			c.mean[t] = c.mean[t] / N;
//		}
//
//		/* il correlatore e' riempito */
//		if ( !c.filled )
//			c.filled = true;
//		/* aggiorno il num. di volte che registro c[N] */
//		c.sweep ++;
//	}
	
	/* crea i "bin" man mano che faccio evolvere il sistema */
//	void oscillatore::update_correlator_binned (void) {
//		/* 'cb.sweep' e' inizializzato a zero */
//		cb.rows = (unsigned int) cb.sweep/cb.nMax + 1;
//		unsigned int k = cb.rows - 1;
//		
//		/* ogni 'cb.nMax' sweep allungo la matrice */
//		if ( cb.sweep % cb.nMax == 0 ) {
//			/* espando 'cb.(**data)' per "righe" */
//			cb.data = (long double **) realloc(cb.data, cb.rows*sizeof(long double *));
//			if ( cb.data == NULL )
//				exit(EXIT_FAILURE);
//
//			/* espando 'cb.(**data)' per "colonne" */
//			cb.data[k] = (long double *) malloc(N*sizeof(long double));
//			if ( cb.data[k] == NULL ) 
//				exit(EXIT_FAILURE);
//
//			/* azzero la matrice */
//			for (unsigned short int j = 0; j < N; j += 1)
//				cb.data[k][j] = (long double) 0;
//		}
//		
//
//		/* aggiorno 'cb.data[k][]' con i correlatori */
//		for (unsigned short int t = 0; t < N; t += 1) {
//			 /* media "ciclica" */
//			for (unsigned short int i = 0; i < N; i += 1)
//				cb.data[k][t] += x[i]*x[(i+t)%N]/N;
//		}
//
//		if ( !cb.initialized )
//			cb.initialized = true;
//		/* incremento numero sweep */
//		cb.sweep++;
//	}
	
	/* normalizza il correlatore (se non lo e') */
//	void oscillatore::normalize_correlator (void) {
//		if ( !c.normalized && c.filled ) {
//			for (unsigned short int t = 0; t < N; t += 1)
//				c.mean[t] = c.mean[t]/c.sweep;
//			c.normalized = true;
//		}
//	}
	
	/* normalizza i bin */
//	void oscillatore::normalize_correlator_binned (void) {
//		if ( !cb.normalized && cb.initialized) {
//			unsigned short int norm = cb.nMax;
//			for (unsigned int k = 0; k < cb.rows; k += 1) {
//				/* calcolo il coefficiente di normalizzazione per l'ultima riga */
//				if ( (k + 1)*cb.nMax > cb.sweep )
//					norm = cb.sweep - k*cb.nMax;
//				/* normalizzo */
//				for (unsigned short int t = 0; t < N; t += 1)
//					cb.data[k][t] = (long double) cb.data[k][t]/norm;
//			}
//			cb.normalized = true;
//		}
//	}


/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: V
 * Description: potenzial energy
 * ------------------------------------------------------------------
 */
long double
oscillatore::V ( long double x ) {
	/* $m \omega^2 x^2 /2 $ */
	return (long double) M * powl( W * x, 2. ) / 2;
//	return 0.; /* particella libera */
} /* -----  end of method oscillatore::V  ----- */

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: diff
 * Description: action difference between two configurations which
 * 				differ only for t-th position
 * ------------------------------------------------------------------
 */
long double
oscillatore::diff ( long double step, unsigned short int t ) {
	/*
	 * XXX le espressioni '( N + t - 1 ) % N' e '( t + 1 ) % N'
	 * servono per assicurare che gli indici del vettore siano
	 * tutti minori di N (e ciclici)
	 */
	long double delta = DELTA * ( 2 * step - 1 );
	long double tmp = 2 * *( x + t ) + delta;
	tmp -= *( x + ( N + t - 1 ) % N ) + *( x + ( t + 1 ) % N );

	return M * delta * tmp / a + a * ( V( *( x + t ) + delta) - V( *( x + t ) ) );
	/*
	 * Lascio il potenziale indicato (non lo svolgo esplicitamente
	 * come ho fatto per l'azione) in modo da poterlo cambiare per 
	 * descrivere altri sistemi con questo algoritmo
	 */
} /* -----  end of method oscillatore::diff  ----- */

	/* riempie il correlatore */
//	void oscillatore::fill_correlator (void) {
//		/* controllo se  e' stato riempito */
//		if ( !c.filled ) {
//			/* se i 'cb.data[][]' esiste, uso quello */
//			if ( cb.initialized ) {
//				/* (eventualmente) normalizzo i bin */
//				normalize_correlator_binned();
//				/* calcolo i correlatori */
//				for (unsigned short int t = 0; t < N; t += 1) {
//					for (unsigned int k = 0; k < cb.rows; k += 1)
//						c.mean[t] += cb.data[k][t];
//					/* normalizzo */
//					c.mean[t] = (long double) c.mean[t]/cb.rows;
//				}
//			}
//			/* altrimenti uso lo storico dei correlatori 'corr.data[][]' */
//			else if ( corr.initialized ) {
//				for (unsigned short int t = 0; t < N; t += 1) {
//					for (unsigned int k = 0; k < corr.rows; k += 1)			
//						c.mean[t] += corr.data[k][t];
//					/* normalizzo */
//					c.mean[t] = c.mean[t]/corr.rows;
//				}
//				/* assegno a 'c.sweep' il numero di sweep */
//				c.sweep = corr.rows;
//			}
//			/* se lo riempio in questa funzione e' anche normalizzato */
//			c.normalized = true;
//		}
//	}
