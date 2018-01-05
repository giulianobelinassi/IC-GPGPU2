#include "shared.h"

#include <cstdio>
#include <cmath>
#include <cublas_v2.h>
#include <thrust/complex.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

extern "C"{

void rigid_body_(
            int* nbe,
            int* n,
            float hest_[],
            float hestdiag[][3][3]
        );

void rigid_body_c_(int* nbe,
		int* n,
		float hest_[],
		float hestdiag[][3][3])
{
	float (*hest)[3*(*nbe)] = (float (*)[3*(*nbe)]) hest_;
	int ma, mb, i, j, ii, jj;

	for (ma = 0; ma < *nbe; ++ma)
	{
		ii = 3*ma;
		for (mb = 0; mb < *n; ++mb)
		{
			jj = 3*mb;
			for (j = 0; j < 3; ++j)
			{
				for (i = 0; i < 3; ++i)
				{
					hestdiag[ma][j][i] -= hest[jj+j][ii+i];
				}
			}

		}
	}

}


__global__ void ghmatece_kernel(
                           int cone[],
                           float cx[],
                           float cy[],
                           float cz[],
                           float cxm[],
                           float cym[],
                           float czm[],
                           float hest[],
                           float rn[][3],
                           float fr,
                           float gi[],
                           float ome[],
                           float c3,
                           float c4,
                           int npg,
                           int n,
                           int nbe,
						   int column_pad,
                           int* ret
                           )
{

	const int ig = threadIdx.y;
	const int jg = threadIdx.x;
	const int ii = blockIdx.y;
	const int jj = blockIdx.x + column_pad;

	const int zgelem_pad = 3*3*npg*npg;
	extern __shared__ float s_[];
	
	int i, j;
	
	float p[4][2], f[4];
	float xj[3][2];
	__shared__ float co[3][4];
	__shared__ float rn_cached[4];
	float p1, p2, p12, g1, g2, sp, sm, rp, rm, det;
	float cxg, cyg, czg, cxp, cyp, czp;
	float j1, j2, j3;
	float r1, r2, r3, r, drn, rd[3];
	float hei;

//Manage the shared memory manually since apparently there is no other way
//to allocate two cubes in dynamic allocated shared memory.
//see https://devblogs.nvidia.com/parallelforall/using-shared-memory-cuda-cc/
#define helem(i, j, k) s_[3*npg*npg*(i) + npg*npg*(j) + (k)]

	int iii, jjj;

	for (iii = 0; iii < 3; ++iii)
		for (jjj = 0; jjj < 3; ++jjj)
			helem(iii, jjj, npg*ig + jg) = 0;

    if (ii != jj)
    {

        if (threadIdx.x < 4 && threadIdx.y == 0)
        {
            co[0][threadIdx.x] = cx[cone[n*threadIdx.x + jj] - 1];
            co[1][threadIdx.x] = cy[cone[n*threadIdx.x + jj] - 1];
            co[2][threadIdx.x] = cz[cone[n*threadIdx.x + jj] - 1];
            rn_cached[threadIdx.x] = rn[jj][threadIdx.x];
        }
        __syncthreads();

        cxp = cxm[ii];
        cyp = cym[ii];
        czp = czm[ii];

        g2 = gi[jg];
        p2 = ome[jg];
        sp = 1 + g2;
        sm = 1 - g2;
        p[0][0] = -0.25f*sm;
        p[1][0] =  0.25f*sm;
        p[2][0] =  0.25f*sp;
        p[3][0] = -0.25f*sp;

        g1 = gi[ig];
        p1 = ome[ig];
        rp = 1 + g1;
        rm = 1 - g1;
        f[0] = 0.25f*rm*sm;
        f[1] = 0.25f*rp*sm;
        f[2] = 0.25f*rp*sp;
        f[3] = 0.25f*rm*sp;
        p[0][1] = -0.25f*rm;
        p[1][1] = -0.25f*rp;
        p[2][1] = 0.25f*rp;
        p[3][1] = 0.25f*rm;

        
       
        for (iii = 0; iii < 2; ++iii)
        {
            for (jjj = 0; jjj < 3; ++jjj)
            {
                xj[jjj][iii] = p[0][iii]*co[jjj][0] + p[1][iii]*co[jjj][1]+ p[2][iii]*co[jjj][2] + p[3][iii]*co[jjj][3];
            }
        }
        

        j1 = xj[1][0]*xj[2][1]-xj[1][1]*xj[2][0];
        j2 = xj[0][1]*xj[2][0]-xj[0][0]*xj[2][1];
        j3 = xj[0][0]*xj[1][1]-xj[0][1]*xj[1][0];

        det = sqrt(j1*j1 + j2*j2 + j3*j3);

        if (det < 1e-5)
        {
            *ret = 1;
            return;
        }


        cxg = 0;
        cyg = 0;
        czg = 0;

        for (iii = 0; iii < 4; ++iii)
        {
            cxg = cxg + co[0][iii]*f[iii];
            cyg = cyg + co[1][iii]*f[iii];
            czg = czg + co[2][iii]*f[iii];
        }

        r1    = cxg - cxp;
        r2    = cyg - cyp;
        r3    = czg - czp;

        r     = sqrt(r1*r1 + r2*r2 + r3*r3);
        drn   = (r1*rn_cached[0] + r2*rn_cached[1] + r3*rn_cached[2])/r;
        rd[0] = r1/r;
        rd[1] = r2/r;
        rd[2] = r3/r;
        
        p12 = p1*p2*det;

        for (j = 0; j < 3; ++j)
        {	for (i = 0; i < 3; ++i)
            {
                
                hei = (c3/(r*r))*(drn*(c4*delta[i][j]+3.0f*rd[i]*rd[j]) +
                    c4*(rd[j]*rn_cached[i]-rd[i]*rn_cached[j]));

                hei = hei*p12;

                helem(j, i, jg*npg + ig) = hei;
            }
        }
        __syncthreads();
    }    
	if (jg < 3 && ig < 3)
	{

		int index = 3*blockIdx.y + (3*nbe)*3*blockIdx.x + ig + (3*nbe)*jg;
		
		hest[index] = thrust::reduce(thrust::seq, &helem(jg, ig, 0), &helem(jg, ig, npg*npg));
	}
}



/*
__global__ void generate_identity(int m, float one_vec[][3])
{
	int i, j;
	const int tid = 3*(blockDim.x*blockIdx.x + threadIdx.x);
	const float const Id[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

	if (tid < m)
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
				one_vec[tid + i][j] = Id[i][j];
	}
}

float* cuda_rigid_body(int nbe, int n, float device_h[])
{
	cudaError_t error;
	cublasStatus_t blaserror;
	float* device_Ids;
	float* device_hdiag;
	const int threads = 32;
	int blocks = (n+threads-1)/threads;
	dim3 threadsPerBlock(32);
	dim3 numBlocks(blocks);
	float one = 1., zero = 0.;	

	error = cudaMalloc(&device_Ids, 3*3*n*sizeof(float));
	cuda_assert(error);
	error = cudaMalloc(&device_hdiag, 3*3*nbe*sizeof(float));
	cuda_assert(error);

	generate_identity<<<numBlocks, threadsPerBlock>>>(n, (float (*)[3]) device_Ids);
	cudaDeviceSynchronize();
	
	cublasHandle_t handle;
	cublasCreate(&handle);

	blaserror = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, nbe, 3, n, &one, device_h, nbe, device_Ids, n, 
			&zero, device_hdiag, nbe);
	cublas_assert(blaserror);

	cudaFree(device_Ids);

	return device_hdiag;
}
*/

__global__ void rigid_body_kernel(int m, int n, float hest_[], float hdiag_[][3][3])
{
	int i, j, k;
	const int tid = blockDim.x*blockIdx.x + threadIdx.x;
	

#define hest(i, j) hest_[(j)*(3*m) + (i)]
#define hdiag(k, i, j) hdiag_[k][j][i]

	if (tid < m)
	{
		for (k = 0; k < n; ++k)
		{
			for (j = 0; j < 3; ++j)
			{
				int jj = 3*k + j;
				for (i = 0; i < 3; ++i)
				{
					int ii = 3*tid + i;
					hdiag(tid, i, j) -= hest(ii, jj);
				}
			}
		}
	}
#undef hest
#undef hdiag
}

float* cuda_rigid_body(int nbe, int n, float device_h[])
{
	cudaError_t error;
	float* device_hdiag;
	const int threads = 128;
	int blocks = (n+threads-1)/threads;
	dim3 threadsPerBlock(threads);
	dim3 numBlocks(blocks);

	error = cudaMalloc(&device_hdiag, 3*3*nbe*sizeof(float));
	cuda_assert(error);

	error = cudaMemset(device_hdiag, 0, 3*3*nbe*sizeof(float));
	cuda_assert(error);

	rigid_body_kernel<<<numBlocks, threadsPerBlock>>>(nbe, n, device_h, (float (*)[3][3]) device_hdiag);

	return device_hdiag;
}


void cuda_ghmatece_(int* nbe,
                    int* npg,
                    int* n,
                    int* np,
                    float* c3,
                    float* c4,
                    float* fr,
                    float hestdiag[][3][3],
                    int* status
                   )
{
	dim3 threadsPerBlock(*npg,*npg);
	int shared_mem_size = 3*3*(*npg)*(*npg)*sizeof(float);
	size_t column_size = (3*(*nbe))*sizeof(float);
	
	cudaError_t error;

	float* device_h;

	int* device_return_status;
	int return_status;
	int width, iterations, i;

    float* hest_ = (float*) malloc(3*(*n)*3*(*nbe)*sizeof(float));
	float (*hest)[3*(*nbe)] = (float (*)[3*(*nbe)]) hest_;

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	width = largest_possible_width(column_size, *n, &iterations);

	error = cudaMalloc(&device_h, (3*(*nbe))*(3*(width))*sizeof(float));
	cuda_assert(error);

	error = cudaMemset(device_return_status, 0, sizeof(int));
	cuda_assert(error);

	for (i = 0; i < iterations; ++i)
	{
		int starting_column = width*i;
		if (starting_column + width > *n)
			width = *n - starting_column;
		dim3 numBlocks(width, *nbe);

		error = cudaMemset(device_h, 0, (3*(*nbe))*(3*(width))*sizeof(float));
		cuda_assert(error);

		cudaDeviceSynchronize();
		ghmatece_kernel<<<numBlocks, threadsPerBlock, shared_mem_size>>>(
							device_cone,
							device_cx,
							device_cy,
							device_cz,
							device_cxm,
							device_cym,
							device_czm,
							device_h,
							(float (*)[3]) device_etas,
							*fr,
							device_gi,
							device_ome,
							*c3,
							*c4,
							*npg,
							*n,
							*nbe,
							starting_column,
							device_return_status
							);
		cudaDeviceSynchronize();

		error = cudaMemcpy(&return_status, device_return_status, sizeof(int), cudaMemcpyDeviceToHost);
		cuda_assert(error);

		if (return_status != 0)
		{
			fputs("Matriz Singular\n", stderr);
		}

//		error = cudaMemcpy(&hest[3*starting_column], device_h, (3*(*nbe))*(3*(width))*sizeof(float), cudaMemcpyDeviceToHost);
//		cuda_assert(error);
	}

	float* device_hdiag = cuda_rigid_body(*nbe, *n, device_h);

	error = cudaFree(device_h);
	cuda_assert(error);
	*status = return_status;
	error = cudaFree(device_return_status);
	cuda_assert(error);

	/*Guarda em shared.cu*/
	device_hestdiag = device_hdiag;

#ifdef TEST_CUDA
	error = cudaMemcpy(hestdiag, device_hdiag, 3*3*(*nbe)*sizeof(float), cudaMemcpyDeviceToHost);
	cuda_assert(error);
#endif
//    rigid_body_(nbe, n, hest_, hestdiag);

}

void cuda_send_gest_data_(int* nbe, float* gestdiag)
{
	cudaError_t error;
	error = cudaMalloc(&device_gestdiag, (*nbe)*3*3*sizeof(float));
	cuda_assert(error);
	error = cudaMemcpy(device_gestdiag, gestdiag, (*nbe)*3*3*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);
}

}
