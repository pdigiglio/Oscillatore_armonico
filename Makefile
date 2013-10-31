# File da compilare
MAIN = oscillatore
INCLUDE = global

# Compilatore c++
CXX = g++

# some cpu-dependent options
MARCH = core2
MASM = intel

# standard language
STD = gnu++11

# Opzioni
CXXFLAGS = -Wall -O2 -Wextra -pedantic -march=$(MARCH) -std=$(STD) \
		   -masm=$(MASM) -mtune=$(MARCH) -fopenmp -lm

# % = $(MAIN); la funzione parte solo 
# se esistono $(MAIN).cpp e 'Makefile'
$(MAIN): %: %.cpp %.cc %.h $(INCLUDE).h Makefile
	@ echo '#INFO'
	@ echo 'Architettura rilevata:\t\t' ` gcc -march=native -Q \
		--help=target | grep --text march | cut -f3 `
	@ echo -e 'Architettura selezionata:\t' $(MARCH)
	@ echo
	$(CXX) $< -o $@ $(CXXFLAGS)
	@ echo
	
# pulisce la directory
clean:
	@ -rm -rf *.d *.o *.tmp $(MAIN)
.PHONY: clean
