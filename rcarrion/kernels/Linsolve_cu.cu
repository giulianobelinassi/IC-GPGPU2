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
//	cublasHandle_t handle;
//	cublasStatus_t stats;
    cudaError_t error;

    thrust::complex<float>* device_zh;
    thrust::complex<float>* device_zfi;
//	cuComplex** device_matrix_pointer;

//    int* device_piv;
//    int* device_status;
	int status;
	int* piv = (int*) malloc((*nn)*sizeof(int));


    thrust::complex<float> one(1.0f, 0.0f);

	error = cudaMalloc(&device_zh, (*nn)*(3*(*n))*sizeof(thrust::complex<float>));
	cuda_assert(error);

    error = cudaMalloc(&device_zfi, (*nn)*sizeof(thrust::complex<float>));
    cuda_assert(error);

//    error = cudaMalloc(&device_piv, (*nn)*sizeof(int));
//   cuda_assert(error);
	
//	error = cudaMalloc(&device_matrix_pointer, sizeof(cuComplex*));
//	cuda_assert(error);

	error = cudaMemcpy(device_zh, zh, (*nn)*(3*(*n))*sizeof(thrust::complex<float>), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_zfi, zfi, (*nn)*sizeof(thrust::complex<float>), cudaMemcpyHostToDevice);
	cuda_assert(error);

	magma_init();
	magma_cgesv_gpu(*nn, 1, (magmaFloatComplex_ptr) device_zh, *nn, piv, (magmaFloatComplex_ptr) device_zfi, *nn, &status);
	magma_finalize();

//	error = cudaMemcpy(device_matrix_pointer, (cuComplex*) &device_zh, sizeof(thrust::complex<float>*), cudaMemcpyHostToDevice);
//	cuda_assert(error);

//    error = cudaMalloc(&device_status, sizeof(int));
//    cuda_assert(error);

//    stats = cublasCreate(&handle);
//	cublas_assert(stats);

//	stats = cublasCgetrfBatched(handle, *nn, (cuComplex**) device_matrix_pointer, *nn, device_piv, device_status, 1);
//    cudaDeviceSynchronize();
//    cublas_assert(stats);

//    error = cudaMemcpy(&status, device_status, sizeof(int), cudaMemcpyDeviceToHost);
//    cuda_assert(error);
//    lu_assert(status);
    
//    stats = cublasCtrsm(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_UNIT, *nn, 1, (cuComplex*) &one, (cuComplex*) device_zh, *nn, (cuComplex*) device_zfi, *nn);
//   cudaDeviceSynchronize();
//    cublas_assert(stats);

//    stats = cublasCtrsm(handle, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, *nn, 1, (cuComplex*) &one, (cuComplex*) device_zh, *nn, (cuComplex*) device_zfi, *nn);
//    cudaDeviceSynchronize();
 //   cublas_assert(stats);

	error = cudaMemcpy(zfi, device_zfi, (*nn)*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
    cuda_assert(error);

    error = cudaFree(device_zh);
    cuda_assert(error);

//    error = cudaFree(device_status);
//    cuda_assert(error);

    error = cudaFree(device_zfi);
    cuda_assert(error);

//    error = cudaFree(device_piv);
//    cuda_assert(error);

//	error = cudaFree(device_matrix_pointer);
//	cuda_assert(error);

	free(piv);
//	cublasDestroy(handle);
//
}
}
