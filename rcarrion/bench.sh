#!/bin/bash

STA_FILES=(ESOLO240E_-5+5.DAT ESOLO960E_-5+5_e10.DAT ESOLO2160E_-5+5.DAT ESOLO4000E_-20+20.DAT ESOLO14400E_-40+40.DAT)
DYN_FILES=(ESOLO240D_-5+5.DAT ESOLO960D_-5+5_e10.DAT ESOLO2160D_-5+5.DAT ESOLO4000D_-20+20.DAT ESOLO14400D_-40+40.DAT)
SOLE_FILES=(SSOLO240E_-5+5.DAT SSOLO960E_-5+5_e10.DAT SSOLO2160E_-5+5.DAT SSOLO4000E_-20+20.DAT SSOLO14400E_-40+40.DAT)
SOLD_FILES=(SSOLO240D_-5+5.DAT SSOLO960D_-5+5_e10.DAT SSOLO2160D_-5+5.DAT SSOLO4000D_-20+20.DAT SSOLO14400D_-40+40.DAT)
MESH_SIZE=(240 960 2160 4000 14400)
COMPILER_PARAMS=(\  gpu gpu \ )
MODE_STR=(gpu gpu_sing cpu)
EXTRA_FLAGS=(\  -f \ )

NUM_WARMUPS=3
NUM_EXECUTIONS=30

for ((i=0; i<3; i++)); do
	compiler_param=${COMPILER_PARAMS[$i]}
	mode_str=${MODE_STR[$i]}
	extra_flag=${EXTRA_FLAGS[$i]}

	make clean
	make ${compiler_param}

	for ((j=0; j<3; ++j)); do
		for thread in {1,8,24,48}; do
			export OMP_NUM_THREADS=${thread}

			for ((k=0;k<5;k++)); do
				mesh=${MESH_SIZE[$k]}
				file_sta=${STA_FILES[$k]}
				file_dyn=${DYN_FILES[$k]}
				sol_sta=${SOLE_FILES[$k]}
				sol_dyn=${SOLD_FILES[$k]}
			
				echo "Executando para $mesh em modo $mode_str com $thread threads"
				echo ""
				for ((l=1;l<=NUM_WARMUPS;l++)); do
					echo "Warmup ${l}: "
					./main $file_sta $file_dyn $sol_sta $sol_dyn $extra_flag
				done
				for ((l=1;l<=NUM_EXECUTIONS;l++)); do
					echo "Execução número $l"
					output_file="results/results_${mode_str}_${mesh}_${thread}_${l}.txt"
					./main $file_sta $file_dyn $sol_sta $sol_dyn $extra_flag > ${output_file}
				done
			done
		done
	done
done

