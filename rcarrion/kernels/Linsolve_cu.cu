#include "shared.h"
#include "../include/magma.h"

#include <thrust/complex.h>
#include <cuda_runtime.h>
#include <cstdio>

extern "C"{

void lu_assert(int status)
{
    if (status != 0)
    {
        fputs("ERRO: Matriz ZH singular\n", stderr);
        exit(1);
    }
}

void cuda_linsolve_(
			int* nn,
			int* n,
			thrust::complex<FREAL> zfi[]
		)
{
    cudaError_t error;

	int status;
	int* piv = (int*) malloc((*nn)*sizeof(int));

    thrust::complex<FREAL> one(1.0f, 0.0f);
    thrust::complex<FREAL> zero(0., 0.);
    
	if (!piv)
	{
		fputs("Erro: Matriz ZH singular", stderr);
		exit(1);
	}

	magma_init();

    if (sizeof(FREAL) == 8)
		magma_zgesv_gpu(*nn, 1, (magmaDoubleComplex_ptr) device_zh, *nn, piv, (magmaDoubleComplex_ptr) device_zfi, *nn, &status);
	else
		magma_cgesv_gpu(*nn, 1, (magmaFloatComplex_ptr) device_zh, *nn, piv, (magmaFloatComplex_ptr) device_zfi, *nn, &status);
	
	magma_finalize();
	lu_assert(status);

	error = cudaMemcpy(zfi, device_zfi, (*nn)*sizeof(thrust::complex<FREAL>), cudaMemcpyDeviceToHost);
    cuda_assert(error);

    error = cudaFree(device_zh);
    cuda_assert(error);

    error = cudaFree(device_zfi);
    cuda_assert(error);

	free(piv);
}
}
