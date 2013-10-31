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
	/* azzero l'azione */
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
 * Description: updates all system positions
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
//	fprintf( stderr, "%u %u: %u\n", c.sweep, c.nB, n );
//	fprintf( stderr, "%Lg\t%Lg\n", c.bin[0][n - 1], c.bin[0][n - 2] );

	/* normalizzo medie ed errori */
	for ( unsigned short int t = 0; t < N; t ++ ) {
		/* normalizzo la media */
		*( c.mean + t ) = *( c.mean + t ) / n;

		/* calcolo l'errore */
		*( c.err + t ) = *( c.err + t ) / n ;
		*( c.err + t ) -= powl( *( c.mean + t ), 2. );
		/* uso la correzione di Bessel per l'errore */
		*( c.err + t ) = sqrtl( *( c.err + t ) / ( n - 1 ) );

		printf( "%hu\t%Lg\t%Lg\n",
				t, *( c.mean + t ) / c.nB, *( c.err + t ) / c.nB );
	}

//	fprintf( stderr, "N. bin: %u\n", n );
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
			*( *( c.bin + t ) + n ) = - *( *( c.bin + t ) + n ) / ( c.nB * ( bns - 1 ) );
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

/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: observables
 * Description: uses clusters to obtain energy density and matrix 
 * 				elemtent value(s)
 * ------------------------------------------------------------------
 */
void
oscillatore::observables ( bool plot = true ) {
	/* variabili temporanee ausiliarie */
	long double ene, tmp;
	/* altre due variabili temporanee */
	long double e_mean, m_mean, e_err, m_err;

	/* numero di bin */
	unsigned int bns = (unsigned) c.sweep / c.nB;

	/* calcolo le energie */
	for ( unsigned short int t = 1; t < 8; t ++ ) {
		/* azzero le varabili temporanee per le medie */
		e_mean = (long double) 0;
		e_err = (long double) 0;

		m_mean = (long double) 0;
		m_err = (long double) 0;

		/* calcolo le medie */
		for ( unsigned int n = 0; n < bns; n ++ ) {
			/* assegno la variabile temporanea */
			ene = *( *( c.bin + t + 1 ) + n ) + *( *( c.bin + t - 1 ) + n );
			ene = acoshl( .5 * ene / *( *( c.bin + t ) + n ) );

			/* aggiorno media ed errore (energia) */
			e_mean += ene;
			e_err += powl( ene, 2. );

			/* assegno la variabile temporanea */
			tmp = expl( - (long double) t * ene );
			tmp = *( *( c.bin + t ) + n ) / ( tmp + expl( - a * N * ene ) / tmp );

			/* aggiorno media ed errore (elemento matrice) */
			m_mean += sqrtl( tmp );
			m_err += tmp;
		}

		/* normalizzo le medie */
		e_mean = e_mean / bns;
		m_mean = m_mean / bns;

		/* calcolo gli errori */
		e_err = e_err / bns - powl( e_mean, 2. );
		e_err = sqrtl( e_err * (bns - 1) );

		m_err = m_err / bns - powl( m_mean, 2. );
		m_err = sqrtl( m_err * (bns - 1) );

		/* stampo a schermo i valori ottenuti in funzione di 't' */
		if( plot ) {
			printf( "%hu\t", t );
			oscillatore::round( e_mean, e_err );
			printf( "\t" );
			oscillatore::round( m_mean, m_err );
			printf( "\n" );
		}
	}
} /* -----  end of method oscillatore::observables  ----- */

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
 * Description: rounds 'val' to the same significant figures  
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

	
/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: V
 * Description: returns potenzial energy
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
