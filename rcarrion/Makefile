MAGMA_INCLUDE = /usr/local/magma/include # Directory containing libmagma header (.h) files
MAGMA_LIB = /usr/local/magma/lib ######### Directory containing libmagma .so file
CUDA_LIB = /usr/local/cuda/lib ########### Directory containing CUDA .so files.

FC      = gfortran-5
NVCC    = nvcc
FCFLAGS = $(EXTRA) -g -cpp -flto -Wall -Wno-unused-dummy-argument -march=native -Ofast -funroll-loops
PARALLEL= -fopenmp
NVFLAGS = $(EXTRA) -I$(MAGMA_INCLUDE) -use_fast_math  -O3  -Xptxas --opt-level=3  -maxrregcount=32 -Xptxas --allow-expensive-optimizations=true 
FLFLAGS = libopenblas.a -fopenmp 
RM      = rm -f

all: main

cpu: main

gpu: MODE = -DUSE_GPU
gpu: FLFLAGS += -lmagma -lcudart -lcublas -lcusparse -lstdc++ -L$(CUDA_LIB) -L$(MAGMA_LIB)
gpu: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Normvec.o Linsolve.o Main.o Ghmatecd_cu.o Ghmatece_cu.o Interec1_cu.o Linsolve_cu.o shared.o 
	$(FC) $(PARALLEL) $(FCFLAGS) -o main $^ $(FLFLAGS)

main: Nonsingd.o Nonsinge.o Sigmaec.o Gauleg.o Sing_de.o Singge.o Solfundif.o Solfune.o Solfund.o Inputece.o Inputecd.o Ghmatece.o Ghmatecd.o Interec.o Outputec.o Normvec.o Linsolve.o Main.o 
	$(FC) $(PARALLEL) $(FCFLAGS) -o $@ $^ $(FLFLAGS)

Main.o: Main.for Inputece.for Ghmatecd.for Ghmatece.for
	$(FC) $(FCFLAGS) $(PARALLEL) $(MODE) -c $<

Linsolve_cu.o: kernels/Linsolve_cu.cu Linsolve.for
	$(NVCC) $(NVFLAGS) -c $<

Ghmatecd_cu.o: kernels/Ghmatecd_cu.cu Sing_de.for Nonsingd.for
	$(NVCC) $(NVFLAGS) -c $<

Ghmatece_cu.o: kernels/Ghmatece_cu.cu kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

Interec1_cu.o: kernels/Interec1_cu.cu kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

shared.o: kernels/shared.cu kernels/shared.h
	$(NVCC) $(NVFLAGS) -c $<

Linsolve.o: Linsolve.for
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Normvec.o: Normvec.for
	$(FC) $(PARALLEL) $(FCFLAGS) -c $<

Inputece.o: Inputece.for
	$(FC) $(FCFLAGS) -c $<

Inputecd.o: Inputecd.for
	$(FC) $(FCFLAGS) -c $<

Ghmatece.o: Ghmatece.for Singge.o Nonsinge.o
	$(FC) $(PARALLEL) $(FCFLAGS) $(MODE) -c $<

Ghmatecd.o: Ghmatecd.for Sing_de.o Nonsingd.o
	$(FC) $(PARALLEL) $(FCFLAGS) $(MODE) -c $<

Interec.o: Interec.for Sigmaec.o Nonsingd.o
	$(FC) $(PARALLEL) $(FCFLAGS) $(MODE) -c $<

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
