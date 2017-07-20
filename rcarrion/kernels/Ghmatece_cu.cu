#include "shared.h"

#include <cstdio>
#include <cmath>
#include <thrust/complex.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

extern "C"{

__global__ void ghmatece_kernel(
                           int cone[],
                           float cx[],
                           float cy[],
                           float cz[],
                           float cxm[],
                           float cym[],
                           float czm[],
                           float hest[],
                           float gest[],
                           float rn[][3],
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
                           int* ret
                           )
{

	const int ig = threadIdx.y;
	const int jg = threadIdx.x;
	const int ii = blockIdx.y;
	const int jj = blockIdx.x;

	int i, j;
	
	float p[4][2], f[4];
	float xj[3][2];
	__shared__ float co[3][4];
	__shared__ float rn_cached[4];
	float p1, p2, p12, g1, g2, sp, sm, rp, rm, det;
	float cxg, cyg, czg, cxp, cyp, czp;
	float j1, j2, j3;
	float r1, r2, r3, r, drn, rd[3];
	float hei, gei;

	__shared__ float helem[3][3][64], gelem[3][3][64];

	int iii, jjj;

	gelem[0][0][npg*ig + jg] = 0;
	gelem[0][1][npg*ig + jg] = 0;
	gelem[0][2][npg*ig + jg] = 0;
	gelem[1][0][npg*ig + jg] = 0;
	gelem[1][1][npg*ig + jg] = 0;
	gelem[1][2][npg*ig + jg] = 0;
	gelem[2][0][npg*ig + jg] = 0;
	gelem[2][1][npg*ig + jg] = 0;
	gelem[2][2][npg*ig + jg] = 0;

	helem[0][0][npg*ig + jg] = 0;
	helem[0][1][npg*ig + jg] = 0;
	helem[0][2][npg*ig + jg] = 0;
	helem[1][0][npg*ig + jg] = 0;
	helem[1][1][npg*ig + jg] = 0;
	helem[1][2][npg*ig + jg] = 0;
	helem[2][0][npg*ig + jg] = 0;
	helem[2][1][npg*ig + jg] = 0;
	helem[2][2][npg*ig + jg] = 0;
	
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
			
			gei = (c1/r)*(c2*delta[i][j]+rd[i]*rd[j]);
			hei = (c3/(r*r))*(drn*(c4*delta[i][j]+3.0f*rd[i]*rd[j]) +
				c4*(rd[j]*rn_cached[i]-rd[i]*rn_cached[j]));

			gei = gei*p12;
			hei = hei*p12;

			gelem[j][i][jg*npg + ig] = gei;
			helem[j][i][jg*npg + ig] = hei;
		}
	}
	__syncthreads();
    
	if (jg < 3 && ig < 3)
	{
		int index = 3*ii + (3*nbe)*3*jj + ig + (3*nbe)*jg;
		
		gest[index] = thrust::reduce(thrust::seq, &gelem[jg][ig][0], &gelem[jg][ig][npg*npg]);
		hest[index] = thrust::reduce(thrust::seq, &helem[jg][ig][0], &helem[jg][ig][npg*npg]);
	}
}

void cuda_ghmatece_(int* nbe,
                    int* npg,
                    int* n,
                    int* np,
                    float* c1,
                    float* c2,
                    float* c3,
                    float* c4,
                    float* fr,
                    float* hest_,
                    float* gest_,
                    int* status
                   )
{
	dim3 threadsPerBlock(*npg,*npg);
	dim3 numBlocks(*n, *nbe);
	cudaError_t error;
    
	float* device_h;
	float* device_g;

	int* device_return_status;
	int return_status;

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	error = cudaMalloc(&device_h, (3*(*nbe))*(3*(*n))*sizeof(float));
	cuda_assert(error);

	error = cudaMalloc(&device_g, (3*(*nbe))*(3*(*n))*sizeof(float));
	cuda_assert(error);

	error = cudaMemset(device_return_status, 0, sizeof(int));
	cuda_assert(error);

	error = cudaMemset(device_h, 0, (3*(*nbe))*(3*(*n))*sizeof(float));
	cuda_assert(error);

	error = cudaMemset(device_g, 0, (3*(*nbe))*(3*(*n))*sizeof(float));
	cuda_assert(error);

	cudaDeviceSynchronize();
	ghmatece_kernel<<<numBlocks, threadsPerBlock>>>(
						device_cone,
						device_cx,
						device_cy,
						device_cz,
						device_cxm,
						device_cym,
						device_czm,
						device_h,
						device_g,
						(float (*)[3]) device_etas,
						*fr,
						device_gi,
						device_ome,
						*c1,
						*c2,
						*c3,
						*c4,
						*npg,
						*n,
						*nbe,
						device_return_status
						);

	cudaDeviceSynchronize();

	error = cudaMemcpy(&return_status, device_return_status, sizeof(int), cudaMemcpyDeviceToHost);
	cuda_assert(error);

	if (return_status != 0)
	{
		fputs("Matriz Singular\n", stderr);
	}

	error = cudaMemcpy(hest_, device_h, (3*(*nbe))*(3*(*n))*sizeof(float), cudaMemcpyDeviceToHost);
	cuda_assert(error);
	error = cudaMemcpy(gest_, device_g, (3*(*nbe))*(3*(*n))*sizeof(float), cudaMemcpyDeviceToHost);
	cuda_assert(error);


	error = cudaFree(device_h);
	cuda_assert(error);
	error = cudaFree(device_g);
	cuda_assert(error);
	*status = return_status;
	error = cudaFree(device_return_status);
	cuda_assert(error);
}
}
