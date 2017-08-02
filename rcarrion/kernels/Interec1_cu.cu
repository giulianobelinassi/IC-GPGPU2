#include "shared.h"

#include <cstdio>
#include <cmath>
#include <thrust/complex.h>
#include <cuda_runtime.h>
//#include "cublas_v2.h"

extern "C"{

/*Reaproveite este kernel pois ele pode ser usado para resolver a parte
 * de deslocamentos de pontos internos de Interec.for*/
__global__ void ghmatecd_kernel(
                           int cone[],
                           float cx[],
                           float cy[],
                           float cz[],
                           float cxm[],
                           float cym[],
                           float czm[],
                           thrust::complex<float> zh[],
                           thrust::complex<float> zg[],
                           float rn[][3],
                           thrust::complex<float> zge,
                           thrust::complex<float> zcs,
                           thrust::complex<float> zcp,
                           float fr,
                           float gi[],
                           float ome[],
                           float c1,
                           float c2,
                           float c3,
                           float c4,
                           int npg,
                           int n,
                           int nbe,
                           int interec,
                           int dim_cone,
                           int* ret
                           );

void cuda_interec1_(int* n,
                    int* nbe,
                    int* npg,
                    int* l,
                    int* np,
					float cxi[],
					float cyi[],
					float czi[],
                    thrust::complex<float>* zge,
                    thrust::complex<float>* zcs,
                    thrust::complex<float>* zcp,
                    float* c1,
                    float* c2,
                    float* c3,
                    float* c4,
                    float* fr,
                    thrust::complex<float> zdfi[],
					thrust::complex<float> zfi[],
                    thrust::complex<float> zdsol[],
                    int* status
                   )
{
	
	dim3 threadsPerBlock(*npg,*npg);
	dim3 numBlocks(*nbe, *l);
	int shared_mem_size = 2*3*3*(*npg)*(*npg)*sizeof(thrust::complex<float>);
	cudaError_t error;
    
	thrust::complex<float> zhelem[3][3];
	thrust::complex<float> zgelem[3][3];

	thrust::complex<float>* device_zh;
	thrust::complex<float>* device_zg;
    thrust::complex<float>* device_zdfi;
    thrust::complex<float>* device_zfi;
    thrust::complex<float>* device_zdsol;

	float* device_cxi;
	float* device_cyi;
	float* device_czi;

	int* device_return_status;
	int return_status;

    thrust::complex<float> one(1., 0.);
    thrust::complex<float> zero(0., 0.);
    thrust::complex<float> minus_one(-1., 0.);

    cublasHandle_t handle;
	cublasStatus_t stats;

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	error = cudaMalloc(&device_zh, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMalloc(&device_zg, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMemset(device_return_status, 0, sizeof(int));
	cuda_assert(error);

	error = cudaMemset(device_zh, 0, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMemset(device_zg, 0, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMalloc(&device_cxi, (*l)*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_cyi, (*l)*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_czi, (*l)*sizeof(float));
	cuda_assert(error);

	error = cudaMemcpy(device_cxi, cxi, (*l)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_cyi, cyi, (*l)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_czi, czi, (*l)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);


	cudaDeviceSynchronize();

	ghmatecd_kernel<<<numBlocks, threadsPerBlock, shared_mem_size>>>(
						device_cone,
						device_cx,
						device_cy,
						device_cz,
						device_cxi,
						device_cyi,
						device_czi,
                        device_zh,
                        device_zg,
						(float (*)[3]) device_etas,
						*zge,
						*zcs,
						*zcp,
						*fr,
						device_gi,
						device_ome,
						*c1,
						*c2,
						*c3,
						*c4,
						*npg,
						*nbe,
						*l,
                        1,
                        *n,
						device_return_status
						);
    
	error = cudaMalloc(&device_zdsol, 3*(*l)*sizeof(thrust::complex<float>));
	cuda_assert(error);

    error = cudaMalloc(&device_zdfi, 3*(*nbe)*sizeof(thrust::complex<float>));
	cuda_assert(error);
	
	error = cudaMalloc(&device_zfi, 3*(*nbe)*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMemcpy(device_zdfi, zdfi, 3*(*nbe)*sizeof(thrust::complex<float>), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_zfi, zfi, 3*(*nbe)*sizeof(thrust::complex<float>), cudaMemcpyHostToDevice);
	cuda_assert(error);
    cudaDeviceSynchronize();

	error = cudaMemcpy(&return_status, device_return_status, sizeof(int), cudaMemcpyDeviceToHost);
	cuda_assert(error);

	if (return_status != 0)
	{
		fputs("Matriz Singular\n", stderr);
	}

    stats = cublasCreate(&handle);
	cublas_assert(stats);
    
	stats = cublasCgemv(handle, CUBLAS_OP_N, 3*(*l), 3*(*nbe), (cuComplex*) &one, (cuComplex*) device_zg, 3*(*l), (cuComplex*) device_zdfi, 1, (cuComplex*) &zero, (cuComplex*) device_zdsol, 1);
    cublas_assert(stats);
	cudaDeviceSynchronize();
    
	stats = cublasCgemv(handle, CUBLAS_OP_N, 3*(*l), 3*(*nbe), (cuComplex*) &(minus_one), (cuComplex*) device_zh, 3*(*l), (cuComplex*) device_zfi, 1, (cuComplex*) &one, (cuComplex*) device_zdsol, 1);
    cublas_assert(stats);
	cudaDeviceSynchronize(); 

	error = cudaMemcpy(zdsol, device_zdsol, 3*(*l)*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
	cuda_assert(error);

    error = cudaFree(device_zh);
	cuda_assert(error);
	error = cudaFree(device_zg);
	*status = return_status;
	error = cudaFree(device_return_status);
	cuda_assert(error);

	error = cudaFree(device_cxi);
	cuda_assert(error);
	error = cudaFree(device_cyi);
	cuda_assert(error);
	error = cudaFree(device_czi);
	cuda_assert(error);
	error = cudaFree(device_zfi);
	cuda_assert(error);
	error = cudaFree(device_zdfi);
	cuda_assert(error);

	error = cudaFree(device_zdsol);
	cuda_assert(error);

	cublasDestroy(handle);
}
}
