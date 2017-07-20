CC      = gcc-5
FC      = gfortran-5
NVCC    = nvcc
FCFLAGS = $(EXTRA) -g -cpp -Wall -Wno-unused-dummy-argument -Ofast -march=native -funroll-loops -flto
PARALLEL= -fopenmp
CFLAGS  = $(EXTRA) -g -funroll-loops -flto -Ofast -Wall -march=native
NVFLAGS = $(EXTRA) -use_fast_math  -O3  -Xptxas --opt-level=3  -maxrregcount=32 -Xptxas --allow-expensive-optimizations=true 
FLFLAGS = libopenblas.a -fopenmp -lcudart -lstdc++ -L/usr/local/cuda-8.0/lib64
RM      = rm -f

all: main

main: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Normvec.o Main.o Ghmatecd_cu.o Ghmatece_cu.o Interec1_cu.o shared.o
	$(FC) $(PARALLEL) $(FCFLAGS) -o $@ $^ $(FLFLAGS)

Main.o: Main.for Inputece.for Ghmatecd.for Ghmatece.for
	$(FC) $(FCFLAGS) -c $<

Ghmatecd_cu.o: kernels/Ghmatecd_cu.cu Sing_de.for Nonsingd.for
	$(NVCC) $(NVFLAGS) -c $<

Ghmatece_cu.o: kernels/Ghmatece_cu.cu kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

Interec1_cu.o: kernels/Interec1_cu.cu kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

shared.o: kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

Normvec.o: Normvec.for
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Inputece.o: Inputece.for
	$(FC) $(FCFLAGS) -c $<

Inputecd.o: Inputecd.for
	$(FC) $(FCFLAGS) -c $<

Ghmatece.o: Ghmatece.for Singge.o Nonsinge.o
	$(FC) $(PARALLEL) $(FCFLAGS) -DGHMATECE_USE_CPU -c $<

Ghmatecd.o: Ghmatecd.for Sing_de.o Nonsingd.o
	$(FC) $(PARALLEL) $(FCFLAGS) -DGHMATECD_USE_CPU -c $<

Interec.o: Interec.for Sigmaec.o Nonsingd.o Interec1_cu.o
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

clean:
	$(RM) main *.o fort.*
