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
			thrust::complex<FREAL> zfi[],
            thrust::complex<FREAL> zdfi[]
		)
{
    cudaError_t error;

	thrust::complex<FREAL>* device_zfi;
    thrust::complex<FREAL>* device_zdfi;

	int status;
	int* piv = (int*) malloc((*nn)*sizeof(int));

    thrust::complex<FREAL> one(1.0f, 0.0f);
    thrust::complex<FREAL> zero(0., 0.);
    
    cublasHandle_t handle;
	cublasStatus_t stats;


	if (!piv)
	{
		fputs("Erro: Matriz ZH singular", stderr);
		exit(1);
	}

	magma_init();

    error = cudaMalloc(&device_zfi, (*nn)*sizeof(thrust::complex<FREAL>));
    cuda_assert(error);

    error = cudaMalloc(&device_zdfi, (*nn)*sizeof(thrust::complex<FREAL>));
    cuda_assert(error);

	error = cudaMemcpy(device_zdfi, zdfi, (*nn)*sizeof(thrust::complex<FREAL>), cudaMemcpyHostToDevice);
	cuda_assert(error);


    stats = cublasCreate(&handle);
	cublas_assert(stats);

    if (sizeof(FREAL) == 8)
        stats = cublasZgemv(handle, CUBLAS_OP_N, (*nn), (*nn), (cuDoubleComplex*) &one, (cuDoubleComplex*) device_zg, (*nn), (cuDoubleComplex*) device_zdfi, 1, (cuDoubleComplex*) &zero, (cuDoubleComplex*) device_zfi, 1); 
    else
        stats = cublasCgemv(handle, CUBLAS_OP_N, (*nn), (*nn), (cuComplex*) &one, (cuComplex*) device_zg, (*nn), (cuComplex*) device_zdfi, 1, (cuComplex*) &zero, (cuComplex*) device_zfi, 1); 
    cublas_assert(stats);
	cublasDestroy(handle);

    error = cudaFree(device_zg);
    cuda_assert(error);
    error = cudaFree(device_zdfi);
    cuda_assert(error);
    
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
