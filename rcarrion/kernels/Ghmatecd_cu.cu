#include "shared.h"

#include <cstdio>
#include <cmath>
#include <thrust/complex.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

extern "C"{

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
						   int column_pad,
                           int* ret
                           )
{

	const int ig = threadIdx.y;
	const int jg = threadIdx.x;
	const int ii = blockIdx.y;
	const int jj = blockIdx.x + column_pad;

	const int gelem_pad = 3*3*npg*npg;
	extern __shared__ thrust::complex<float> s_[];
	
	int i, j;
	
	const float pi  = 3.141592654;
	float p[4][2], f[4];
	float xj[3][2];
	__shared__ float co[3][4];
	__shared__ float rn_cached[4];
	float g1, g2, p1, p2, p12, sp, sm, rp, rm, det;
	float cxg, cyg, czg, cxp, cyp, czp;
	float j1, j2, j3;
	float r1, r2, r3, r, drn, rd[3];
	thrust::complex<float> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                    zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	thrust::complex<float> zhi, zgi;

	int iii, jjj;

#define zhelem(i, j, k) s_[3*npg*npg*(i) + npg*npg*(j) + (k)]
#define zgelem(i, j, k) (s_ + gelem_pad)[3*npg*npg*(i) + npg*npg*(j) + (k)]

	zgelem(0, 0, npg*ig + jg) = 0;
	zgelem(0, 1, npg*ig + jg) = 0;
	zgelem(0, 2, npg*ig + jg) = 0;
	zgelem(1, 0, npg*ig + jg) = 0;
	zgelem(1, 1, npg*ig + jg) = 0;
	zgelem(1, 2, npg*ig + jg) = 0;
	zgelem(2, 0, npg*ig + jg) = 0;
	zgelem(2, 1, npg*ig + jg) = 0;
	zgelem(2, 2, npg*ig + jg) = 0;

	zhelem(0, 0, npg*ig + jg) = 0;
	zhelem(0, 1, npg*ig + jg) = 0;
	zhelem(0, 2, npg*ig + jg) = 0;
	zhelem(1, 0, npg*ig + jg) = 0;
	zhelem(1, 1, npg*ig + jg) = 0;
	zhelem(1, 2, npg*ig + jg) = 0;
	zhelem(2, 0, npg*ig + jg) = 0;
	zhelem(2, 1, npg*ig + jg) = 0;
	zhelem(2, 2, npg*ig + jg) = 0;
	
	if (threadIdx.x < 4 && threadIdx.y == 0)
	{
		co[0][threadIdx.x] = cx[cone[dim_cone*threadIdx.x + jj] - 1];
		co[1][threadIdx.x] = cy[cone[dim_cone*threadIdx.x + jj] - 1];
		co[2][threadIdx.x] = cz[cone[dim_cone*threadIdx.x + jj] - 1];
		//Note que a dimensão coluna de rn é 3, mas estamos acessando o elemento
		//na posição 4. Isto pode levar a um segfault, entretanto consegue-se
		//uma melhora de ~100ms no kernel se fizermos esta alteração.
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
	

	zwi = thrust::complex<float>(0, fr);
	
	zc0 = 1.f/(4.f*(pi)*(zge));
	zc1 = ((zcp)/(zcs))*((zcp)/(zcs));
	zc2 = ((zcs)/(zcp))*((zcs)/(zcp));
	zkp = -zwi/(zcp);
	zks = -zwi/(zcs);
	zzp = zkp*r;
	zzs = zks*r;
	zezp= exp(zzp);
	zezs= exp(zzs);
	zp2 = zzp*zzp;
	zs2 = zzs*zzs;

	zfhi    = (1.f + 1.f/zs2 - 1.f/zzs)*zezs/r - zc2*(1.f/zp2 - 1.f/zzp)*zezp/r;
	zcappa  = (1.f + 3.f/zs2 - 3.f/zzs)*zezs/r - zc2*(1.f + 3.f/zp2 - 3.f/zzp)*zezp/r;
	zfhidr  = (-2.f+ zzs + 3.f/zzs - 3.f/zs2)*zezs/(r*r) - zc2*(-1.f + 3.f/zzp - 3.f/zp2)*zezp/(r*r);
	zcappadr= (zzs - 4.f + 9.f/zzs - 9.f/zs2)*zezs/(r*r) - zc2*(zzp - 4.f + 9.f/zzp - 9.f/zp2)*zezp/(r*r);

	zaa = zfhidr-zcappa/r;
	zbb = 4.f*zcappa/r - 2.f*zcappadr;
	zcc = (zc1-2.f)*(zaa + 0.5f*zbb-3.0f*zcappa/r)-2.0f*zcappa/r;

	p12 = p1*p2*det;
	
    if (ii != jj)
    {
        for (j = 0; j < 3; ++j)
        {	for (i = 0; i < 3; ++i)
            {
                zgi = (zc0*(zfhi*delta[j][i] - zcappa*rd[j]*rd[i]));
                

                zhi = (1.0f/(4.0f*pi))*((zaa*(drn*delta[j][i] + 
                                    rd[j]*rn_cached[i])) + rd[i]*rd[j]*drn*zbb + 
                            rd[i]*rn_cached[j]*zcc);
            
                if (ii == jj && !interec)
                {
                    zgi = zgi - (c1/r)*(c2*delta[j][i] + rd[i]*rd[j]);
                    zhi = zhi - (c3/(r*r))*(drn*(c4*delta[j][i] + 3.0f*rd[i]*rd[j]) + c4*(rd[j]*rn_cached[i] - rd[i]*rn_cached[j]));
                }
                
                zgi = zgi*p12;
                zhi = zhi*p12;

                zgelem(j, i, jg*npg + ig) = zgi;
                zhelem(j, i, jg*npg + ig) = zhi;
            }
        }
    }
	__syncthreads();
	
	if (jg < 3 && ig < 3)
	{
		int index = 3*blockIdx.y + (3*nbe)*3*blockIdx.x + ig + (3*nbe)*jg;

		zg[index] = thrust::reduce(thrust::seq, &zgelem(jg, ig, 0), &zgelem(jg, ig, npg*npg));
		zh[index] = thrust::reduce(thrust::seq, &zhelem(jg, ig, 0), &zhelem(jg, ig, npg*npg));
	}
}


void cuda_ghmatecd_(int* nbe,
                    int* npg,
                    int* n,
                    int* np,
                    thrust::complex<float>* zge,
                    thrust::complex<float>* zcs,
                    thrust::complex<float>* zcp,
                    float* c1,
                    float* c2,
                    float* c3,
                    float* c4,
                    float* fr,
                    float* zhest_,
                    float* zgest_,
                    thrust::complex<float>* zgp_,
                    thrust::complex<float>* zhp_,
                    int* status
                   )
{
	
	size_t column_size = 2*(3*(*nbe))*sizeof(thrust::complex<float>);
	
	int shared_mem_size = 2*3*3*(*npg)*(*npg)*sizeof(thrust::complex<float>);
	cudaError_t error;
	
	thrust::complex<float>* device_zh;
	thrust::complex<float>* device_zg;

	int* device_return_status;
	int return_status;

	/*Cast os parâmetros de volta para o tipo original*/
	float (*zhest)[3*(*nbe)] = (float (*)[3*(*nbe)]) zhest_;
	float (*zgest)[3*(*nbe)] = (float (*)[3*(*nbe)]) zgest_;
	thrust::complex<float> (*zgp)[3*(*nbe)] = (thrust::complex<float> (*)[3*(*nbe)]) zgp_;
	thrust::complex<float> (*zhp)[3*(*nbe)] = (thrust::complex<float> (*)[3*(*nbe)]) zhp_;

	int i, ii, iterations, width;
	dim3 threadsPerBlock(*npg,*npg);

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	width = largest_possible_width(column_size, *n, &iterations);

	error = cudaMalloc(&device_zh, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMalloc(&device_zg, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<float>));
	cuda_assert(error);

	for (i = 0; i < iterations; ++i)
	{
		int starting_column = width*i;
		if (starting_column + width > *n)
			width = *n - starting_column;
		dim3 numBlocks(width, *nbe);


		error = cudaMemset(device_return_status, 0, sizeof(int));
		cuda_assert(error);

		error = cudaMemset(device_zh, 0, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<float>));
		cuda_assert(error);

		error = cudaMemset(device_zg, 0, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<float>));
		cuda_assert(error);

		cudaDeviceSynchronize();
		ghmatecd_kernel<<<numBlocks, threadsPerBlock, shared_mem_size>>>(
							device_cone,
							device_cx,
							device_cy,
							device_cz,
							device_cxm,
							device_cym,
							device_czm,
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
							*n,
							*nbe,
							0,
							*n,
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

		error = cudaMemcpy(&zhp[3*starting_column], device_zh, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
		cuda_assert(error);
		error = cudaMemcpy(&zgp[3*starting_column], device_zg, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
		cuda_assert(error);

	}

	error = cudaFree(device_zh);
	cuda_assert(error);
	error = cudaFree(device_zg);
	cuda_assert(error);
	*status = return_status;
    error = cudaFree(device_return_status);
    cuda_assert(error);
}
}
