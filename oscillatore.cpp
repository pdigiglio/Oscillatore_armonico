/* Contiene l'implementazione della classe 'oscillatore' */
#include "oscillatore.cc"

/* per il seme dei numeri casuali */
#include <time.h>

int main (void) {
	/* inizializzo generatore num. casuali */
	srand( time( NULL ) );
		
	oscillatore sistema;

	/* termalizzazione */
	for ( unsigned short int i = 0; i < 100; i ++ )
		sistema.sweep();

	/* misure */
	for ( unsigned int i = 0; i < 100; i ++ ) {
		sistema.sweep();
		sistema.correlator();
		sistema.autocorrelator(0);
	}

	sistema.plot_autocorrelator();

//	sistema.print_c();
//	sistema.Cluster();
//	sistema.observables();

	/*
	sistema.Cluster();
//	printf( " > ENERGIE:\n" );
	fprintf( stdout, "#t\tene\t\terr\t\tmat\t\terr\n" );
	sistema.observables();
//	sistema.print_c();

//	printf( "\n > ELEMENTO MATRICE:\n" );
//	sistema.matrix();
	*/
	exit(EXIT_SUCCESS);
}
