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
			thrust::complex<float> zh[],
			thrust::complex<float> zfi[]
		)
{
    cudaError_t error;

    magmaFloatComplex_ptr device_zh;
    magmaFloatComplex_ptr device_zfi;
	int status;
	int* piv = (int*) malloc((*nn)*sizeof(int));

    thrust::complex<float> one(1.0f, 0.0f);

	if (!piv)
	{
		fputs("Erro: Matriz ZH singular", stderr);
		exit(1);
	}

	magma_init();

	error = cudaMalloc(&device_zh, (*nn)*(*nn)*sizeof(thrust::complex<float>));
	cuda_assert(error);

    error = cudaMalloc(&device_zfi, (*nn)*sizeof(thrust::complex<float>));
    cuda_assert(error);

	error = cudaMemcpy(device_zh, zh, (*nn)*(*nn)*sizeof(thrust::complex<float>), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_zfi, zfi, (*nn)*sizeof(thrust::complex<float>), cudaMemcpyHostToDevice);
	cuda_assert(error);

	magma_cgesv_gpu(*nn, 1, device_zh, *nn, piv, device_zfi, *nn, &status);
	magma_finalize();
	lu_assert(status);

	error = cudaMemcpy(zfi, device_zfi, (*nn)*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
    cuda_assert(error);

    error = cudaFree(device_zh);
    cuda_assert(error);

    error = cudaFree(device_zfi);
    cuda_assert(error);

	free(piv);
}
}
