#include "shared.h"

#include <cstdio>
#include <cmath>
#include <thrust/complex.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

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
                           int* ret
                           );

void cuda_interec1_(int* nbe,
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
                    thrust::complex<float>* zgp_,
                    thrust::complex<float>* zhp_,
                    thrust::complex<float> zdfi[],
					thrust::complex<float> zfi[],
                    int* status
                   )
{
	dim3 threadsPerBlock(*npg,*npg);
	dim3 numBlocks(*nbe, *l);
	cudaError_t error;
    
	thrust::complex<float> zhelem[3][3];
	thrust::complex<float> zgelem[3][3];

	thrust::complex<float>* device_zh;
	thrust::complex<float>* device_zg;
	float* device_cxi;
	float* device_cyi;
	float* device_czi;

	int* device_return_status;
	int return_status;

	thrust::complex<float> (*zgp)[3*(*l)] = (thrust::complex<float> (*)[3*(*l)]) zgp_;
	thrust::complex<float> (*zhp)[3*(*l)] = (thrust::complex<float> (*)[3*(*l)]) zhp_;

	int i, ii;

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

	ghmatecd_kernel<<<numBlocks, threadsPerBlock>>>(
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
						device_return_status
						);
	cudaDeviceSynchronize();

	error = cudaMemcpy(&return_status, device_return_status, sizeof(int), cudaMemcpyDeviceToHost);
	cuda_assert(error);

	if (return_status != 0)
	{
		fputs("Matriz Singular\n", stderr);
	}

	error = cudaMemcpy(zhp_, device_zh, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
	cuda_assert(error);
	error = cudaMemcpy(zgp_, device_zg, (3*(*nbe))*(3*(*l))*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
	cuda_assert(error);
/*
	for (int j = 0; j < 3*(*nbe); ++j)
	{
		for (int i = 0; i < 3*(*l); ++i)
		{
			printf("%d %d %f %f\n", i, j, zgp[j][i].real(), zgp[j][i].imag());
		}
	}
*/
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
}
}
