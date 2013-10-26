/* Contiene l'implementazione della classe 'oscillatore' */
#include "oscillatore.cc"

/* per il seme dei numeri casuali */
#include <time.h>

int main (void) {
	/* inizializzo generatore num. casuali */
	srand( time( NULL ) );
		
	oscillatore sistema;
//	sistema.plot_state();
	for (unsigned int i = 0; i < 1000; i += 1) {
//		printf("%u\t%Lf\n", i, sistema.get_action());
		sistema.sweep();
		if ( i >= 100 )
			sistema.update_correlator_binned();
	}
	
	sistema.normalize_correlator_binned();
//	sistema.plot_correlator();
//	sistema.normalize_correlator_binned();
	
//	sistema.autocorrelator(1);
//	sistema.plot_autocorrelator();

//	sistema.set_nMax(150);
	sistema.correlator_bin();

	sistema.create_clusters();
//	sistema.correlator_errors(true);
	
	/* se passo 'true' come argomento stampa un grafico */
	printf("Energie:\n");
	sistema.energy(true);
	printf("Elemento matrice:\n");
	sistema.ananas(true);
	
	exit(EXIT_SUCCESS);
}
