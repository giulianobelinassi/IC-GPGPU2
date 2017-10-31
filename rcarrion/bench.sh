#!/bin/bash

make clean
make gpu EXTRA="-DTEST_CUDA" -j 4

#warmup step
for i in {1..3}; do
./main
done

for i in {1..1}; do
	echo "Executando em paralelo com 240: $i"
	./main > "results/parallel_240_$i.txt"

done

for i in {1..1}; do
	echo "Executando em paralelo com 960: $i"
	./main ESOLO960E_-5+5_e10.DAT ESOLO960D_-5+5_e10.DAT SSOLO960E_-5+5_e10.DAT SSOLO960D_-5+5_e10.DAT > "results/parallel_960_$i.txt"
done
for i in {1..1}; do
	echo "Executando em paralelo com 2160: $i"
	./main ESOLO2160E_-5+5.DAT ESOLO2160D_-5+5.DAT SSOLO2160E_-5+5.DAT SSOLO2160D_-5+5.DAT > "results/parallel_2160_$i.txt"
done
for i in {1..1}; do
	echo "Executando em paralelo com 4000: $i"
	./main ESOLO4000E_-20+20.DAT ESOLO4000D_-20+20.DAT SSOLO4000E_-20+20.DAT SSOLO4000D_-20+20.DAT > "results/parallel_4000_$i.txt"
done

make clean
make cpu

OMP_NUM_THREADS=1
#warmup
for i in {1..1}; do
	./main
done

for i in {1..1}; do
	echo "Executando em série com 240: $i"
	./main > "results/sequential_240_$i.txt"
done

for i in {1..1}; do
	echo "Executando em série com 960: $i"
	./main ESOLO960E_-5+5_e10.DAT ESOLO960D_-5+5_e10.DAT SSOLO960E_-5+5_e10.DAT SSOLO960D_-5+5_e10.DAT > "results/sequential_960_$i.txt"
done

for i in {1..1}; do
	echo "Executando em série com 2160: $i"
	./main ESOLO2160E_-5+5.DAT ESOLO2160D_-5+5.DAT SSOLO2160E_-5+5.DAT SSOLO2160D_-5+5.DAT > "results/sequential_2160_$i.txt"
done

for i in {1..1}; do
	echo "Executando em série com 4000: $i"
	./main ESOLO4000E_-20+20.DAT ESOLO4000D_-20+20.DAT SSOLO4000E_-20+20.DAT SSOLO4000D_-20+20.DAT > "results/sequential_4000_$i.txt"
done


