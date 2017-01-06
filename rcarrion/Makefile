CC      = gcc
FC      = gfortran
FCFLAGS = -g -std=legacy -flto -Wall -Wno-unused-dummy-argument -march=native -Ofast
PARALLEL= -fopenmp  
CFLAGS  = -O3 -march=native
FLFLAGS = libopenblas.a
RM      = rm -f

all: main main_small

main: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Main.o unl_stack.o
	$(FC) $(PARALLEL) $(FCFLAGS) -o $@ $^ $(FLFLAGS)

main_small: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Main_small.o unl_stack.o
	$(FC) $(PARALLEL) $(FCFLAGS) -o $@ $^ $(FLFLAGS)

Main.o: Main.for Inputece.for
	$(FC) $(FCFLAGS) -c $<

Main_small.o: Main_small.for Inputece.for
	$(FC) $(FCFLAGS) -c $<

Inputece.o: Inputece.for
	$(FC) $(FCFLAGS) -c $<

Inputecd.o: Inputecd.for
	$(FC) $(FCFLAGS) -c $<

Ghmatece.o: Ghmatece.for Singge.o Nonsinge.o
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Ghmatecd.o: Ghmatecd.for Sing_de.o Nonsingd.o
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Interec.o: Interec.for Sigmaec.o Nonsingd.o
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Outputec.o: Outputec.for
	$(FC) $(FCFLAGS) -c $<

Sigmaec.o: Sigmaec.for Gauleg.o
	$(FC) $(FCFLAGS) -c $<

Sing_de.o: Sing_de.for Gauleg.o
	$(FC) $(FCFLAGS) -c $<

Singge.o: Singge.for Nonsinge.o
	$(FC) $(FCFLAGS) -c $<

Solfund.o: Solfund.for
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Solfundif.o: Solfundif.for
	$(FC) $(FCFLAGS) -c $<

Solfune.o: Solfune.for
	$(FC) $(FCFLAGS) -c $<

Gauleg.o: Gauleg.for
	$(FC) $(FCFLAGS) -c $<

Nonsingd.o: Nonsingd.for
	$(FC) $(FCFLAGS) -c $<

Nonsinge.o: Nonsinge.for
	$(FC) $(FCFLAGS) -c $<

unl_stack.o: unl_stack.c
	$(CC) $(CFLAGS) -c $<
	

clean:
	$(RM) main main_small *.o fort.*
