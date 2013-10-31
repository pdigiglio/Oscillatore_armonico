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
	for ( unsigned int i = 0; i < 10000; i ++ ) {
		sistema.sweep();
		sistema.correlator();
//		sistema.autocorrelator(0);
	}

	sistema.Cluster();
	printf( "> ENERGIE:\n" );
	sistema.energy();
//	sistema.print_c();

	printf( "> ELEMENTO MATRICE:\n" );
	sistema.matrix();
	
	exit(EXIT_SUCCESS);
}
