FC      = gfortran
FCFLAGS = -g -std=legacy #-Wall
FLFLAGS = -llapack
RM      = rm -f

all: main

main: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Main.o
	$(FC) $(FCFLAGS) -o $@ $^ $(FLFLAGS)

Main.o: Main.for Inputece.for
	$(FC) $(FCFLAGS) -c $<

Inputece.o: Inputece.for
	$(FC) $(FCFLAGS) -c $<

Inputecd.o: Inputecd.for
	$(FC) $(FCFLAGS) -c $<

Ghmatece.o: Ghmatece.for Singge.o Nonsinge.o
	$(FC) $(FCFLAGS) -c $<

Ghmatecd.o: Ghmatecd.for Sing_de.o Nonsingd.o
	$(FC) $(FCFLAGS) -c $<

Interec.o: Interec.for Sigmaec.o Nonsingd.o
	$(FC) $(FCFLAGS) -c $<

Outputec.o: Outputec.for
	$(FC) $(FCFLAGS) -c $<

Sigmaec.o: Sigmaec.for Gauleg.o
	$(FC) $(FCFLAGS) -c $<

Sing_de.o: Sing_de.for Gauleg.o
	$(FC) $(FCFLAGS) -c $<

Singge.o: Singge.for Nonsinge.o
	$(FC) $(FCFLAGS) -c $<

Solfund.o: Solfund.for
	$(FC) $(FCFLAGS) -c $<

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



clean:
	$(RM) main *.o fort.*
