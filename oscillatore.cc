/* C standard math library */
#include <math.h>

/* Class header file */
#include "oscillatore.h"

/* colors */
#define ANSI_RED     "\x1b[31m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_YELLOW  "\x1b[33m"
#define ANSI_BLUE    "\x1b[34m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_CYAN    "\x1b[36m"
#define ANSI_RESET   "\x1b[0m"
#define ANSI_BOLD	"\033[1m"

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
		*( bin + t ) = (long double *) malloc( sizeof(long double) );
		/* alloco la memoria temporanea per gli autocorrelatori */
		*( corr + t ) = (long double *) malloc( sizeof(long double) );	
		/* prova di output colorato */
		if ( *( bin + t ) == NULL || *( corr + t ) == NULL ) {
			fprintf ( stderr,
					"[" ANSI_BOLD ANSI_MAGENTA "%s" ANSI_RESET "]"
					" dynamic memory allocation failed.\n",
					(char *) __PRETTY_FUNCTION__ );
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
		*( bin + t ) = (long double *) malloc( sizeof(long double) );
		if ( *( bin + t ) == NULL ) {
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

	void oscillatore::plot_autocorrelator( void ) {
			mean = mean * mean;
			/* lo stampo normalizzato a 'ac[0]' */
			for ( unsigned short int j = 0; j < 30; j ++ )
				printf("%u\t%Lf\n", j, ( *( ac + j ) - mean ) * ( 1. + (long double) j / ( c.sweep - j) ) / ( *ac - mean ) );
	}

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
	bs = size;	
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
	unsigned int b = (unsigned) c.sweep / bs;

//	printf( "Sweep: %u\n", c.sweep );

	/* calcolo i correlatori */
	for ( unsigned short int t = 0; t < N; t ++ ) {
		/* azzero la variabile ausiliaria */
		*( *( corr + t) + c.sweep % count ) = (long double) 0;
		/* calcolo il t-esimo correlatore */
		for ( unsigned short int i = 0; i < N; i ++ )
			*( *( corr + t) + c.sweep % count ) += *( x + i ) * *( x + ( i + t ) % N );

		/* aggiorno il bin b-esimo del correlatore t-esimo */
//		*( *( corr + t ) + c.sweep % count ) = *( *( corr + t ) + c.sweep % count ) / N;
		*( *( bin + t ) + b ) += *( *( corr + t ) + c.sweep % count );
	}
	
	/*
	 * Aggiorno il numero di sweep e controllo se devo allungare la
	 * matrice dei bin. Questa matrice viene allungata quando tutti
	 * gli elementi sono stati sistemati in un bin
	 */
	if ( !( ++ c.sweep % bs ) ) {
		/* ricalcolo il valore di 'b' */
		b = (unsigned) c.sweep / bs;
		/* alloco la memoria per i nuovi bin */
		for ( unsigned short int t = 0; t < N; t ++ ) {
			*( bin + t ) = (long double *) realloc( *( bin + t), ( b + 1 ) * sizeof(long double) );
			/* controllo che non ci siano problemi di memoria */
			if ( bin + t == NULL ) {
				fprintf ( stderr, "[correlator] Dynamic memory reallocation failed\n" );
				exit(EXIT_FAILURE);
			}
			
			/* azzero le variabili (bin successivi) */
			*( *( bin + t ) + b ) = (long double) 0;

			/* aggiorno la media t-esima */
			*( c.mean + t ) += *( *( bin + t ) + b - 1 );
			/* aggiorno l'errore */
			*( c.err + t ) += powl( *( *( bin + t ) + b - 1 ), 2. );
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
	unsigned int n = (unsigned) c.sweep / bs;
//	fprintf( stderr, "%u %u: %u\n", c.sweep, bs, n );
//	fprintf( stderr, "%Lg\t%Lg\n", bin[0][n - 1], bin[0][n - 2] );

	/* normalizzo medie ed errori */
	for ( unsigned short int t = 0; t < N; t ++ ) {
		/* normalizzo la media */
		*( c.mean + t ) = *( c.mean + t ) / n;

		/* calcolo l'errore */
		*( c.err + t ) = *( c.err + t ) / n ;
		*( c.err + t ) -= powl( *( c.mean + t ), 2. );
		/* uso la correzione di Bessel per l'errore */
		*( c.err + t ) = sqrtl( *( c.err + t ) / ( n - 1 ) );

		printf( "%hu\t", t );
		oscillatore::round( *( c.mean + t ) / bs, *( c.err + t ) / bs );
		printf( "\n" );
	}

//	fprintf( stderr, "N. bin: %u\n", n );
} /* -----  end of method oscillatore::print_c  ----- */


/*
 * ------------------------------------------------------------------
 *       Class: oscillatore
 *      Method: cluster
 * Description: creates non-normalized clusters from correlators bin
 * ------------------------------------------------------------------
 */
void
oscillatore::Cluster ( void ) {
	/* numero di bin */
	unsigned int bns = (unsigned) c.sweep / bs;

	/* assegno alle variabili contenenti i bin i cluster rispettivi */
	for ( unsigned int n = 0; n < bns; n ++ )
		for ( unsigned short int t = tMin - 1; t <= tMax + 1; t ++ ) {
			*( *( bin + t ) + n ) -= *( c.mean + t );
			*( *( bin + t ) + n ) = - *( *( bin + t ) + n ); // / ( N * bs * ( bns - 1 ) );
		}
} /* -----  end of method oscillatore::cluster  ----- */

	
	/* 
	 * calcola l'autocorrelatore
	 * stop = num. a cui mi fermo (es. ne calcolo ~ 100)
	 * t = tempo reticolo
	 */
	void oscillatore::autocorrelator ( unsigned short int t = 0 ) {
//		/* FIXME */
//		for ( unsigned short int k = 0; k <= MIN( c.sweep, count ); k ++ ) {
//			*( ac + k ) += *( *( corr + t) + c.sweep % count ) * *( *( corr + t ) + ( count + c.sweep - k ) % count );
////			printf( "%hu/%hu:\t%Lg\n", k, MIN( c.sweep - 1, 30 ), *( ac + k ) );
//		}
//
//		mean += *( *(corr + t) );
//
//
//		count ++;
//		for ( unsigned short j = 0; j < N; j ++ ) {
//
//			*( corr + j ) = (long double *) realloc ( *( corr + j ), count * sizeof(long double) );
//			if ( *( corr + j ) == NULL ) {
//				fprintf ( stderr, "\ndynamic memory reallocation failed\n" );
//				exit (EXIT_FAILURE);
//			}
//
//		}
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
	/* altre variabili temporanee */
	long double e_mean, m_mean, e_err, m_err;

	/* numero di bin */
	unsigned int bns = (unsigned) c.sweep / bs;

	/* calcolo le energie */
	for ( unsigned short int t = tMin; t <= tMax; t ++ ) {
		/* azzero le varabili temporanee per le medie */
		e_mean = (long double) 0;
		e_err = (long double) 0;

		m_mean = (long double) 0;
		m_err = (long double) 0;

		/* calcolo le medie */
		for ( unsigned int n = 0; n < bns; n ++ ) {
			/* TODO implementa la creazione dei cluster qui dentro */
//			*( *( bin + t ) + n ) -= *( c.mean + t );
//			*( *( bin + t ) + n ) = - *( *( bin + t ) + n );
			/* assegno la variabile temporanea */
			ene = *( *( bin + t + 1 ) + n ) + *( *( bin + ( N + t - 1 ) % N ) + n );
			ene = acoshl( .5 * ene / *( *( bin + t ) + n ) );

			/* aggiorno media ed errore (energia) */
			e_mean += ene;
			e_err += powl( ene, 2. );

			/* assegno la variabile temporanea */
			tmp = expl( - (long double) t * ene );
			tmp = *( *( bin + t ) + n ) / ( tmp + expl( - a * N * ene ) / tmp );

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

		/* 
		 * :REMARK:03/11/2013 23:17:55:: Sistemo la normalizzazione
		 * che ho trascurato nella creazione dei cluster
		 */
		m_err = m_err / bns - powl( m_mean, 2. );
		m_err = sqrtl( m_err / ( N * bs ) );

		/* stampo a schermo i valori ottenuti in funzione di 't' */
		if( plot /* && t > 0 && t < 8 */ ) {
			printf( "%hu\t", t );
			oscillatore::round( e_mean, e_err );
			printf( "\t" );	
			/* 
			 * :REMARK:03/11/2013 23:17:55:: Sistemo la normalizzazione
			 * che ho trascurato nella creazione dei cluster
			 */
			oscillatore::round( m_mean / sqrtl( N * bs * ( bns - 1 ) ), m_err );
			printf( "\n" );
		}
	}
} /* -----  end of method oscillatore::observables  ----- */

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
				floorl( val / powl(10., exp) + 0.5) * powl(10., exp),
				floorl( err / powl(10., exp) + 0.5) * powl(10., exp)
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
	 * XXX le espressioni '( N + t - 1 ) % N' e '( t + 1 ) % N' servo-
	 * no per assicurare che gli indici del vettore siano tutti minori
	 * di N (e ciclici)
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
