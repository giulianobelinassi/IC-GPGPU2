#include "shared.h"

#include <cstdio>
#include <cmath>
#include <thrust/complex.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

extern "C"{

__global__ void ghmatecd_kernel(
                           int cone[],
                           FREAL cx[],
                           FREAL cy[],
                           FREAL cz[],
                           FREAL cxm[],
                           FREAL cym[],
                           FREAL czm[],
                           thrust::complex<FREAL> zh[],
                           thrust::complex<FREAL> zg[],
                           FREAL rn[][3],
                           thrust::complex<FREAL> zge,
                           thrust::complex<FREAL> zcs,
                           thrust::complex<FREAL> zcp,
                           FREAL fr,
                           FREAL gi[],
                           FREAL ome[],
                           FREAL c1,
                           FREAL c2,
                           FREAL c3,
                           FREAL c4,
//                           int npg,
                           int n,
                           int nbe,
                           int dim_cone,
						   int column_pad,
                           int fast_singular,
                           int* ret
                           )
{

	const int ig = threadIdx.y;
	const int jg = threadIdx.x;
	const int ii = blockIdx.y;
	const int jj = blockIdx.x + column_pad;

	const int npg = blockDim.y;

	const int tid = npg*threadIdx.y + threadIdx.x;
	const int lane = tid % 32;
	const int warp = tid / 32;
	const int num_warps = (npg*npg + 31)/32;

	const int lane_x = lane % 3;
	const int lane_y = lane / 3;

	const int gelem_pad = 3*3*num_warps;
	extern __shared__ thrust::complex<FREAL> s_[];
	
	int i, j;
	
	const FREAL pi  = 3.141592654;
	FREAL p[4][2], f[4];
	FREAL xj[3][2];
	__shared__ FREAL co[3][4];
	__shared__ FREAL rn_cached[4];
	FREAL g1, g2, p1, p2, p12, sp, sm, rp, rm, det;
	FREAL cxg, cyg, czg, cxp, cyp, czp;
	FREAL j1, j2, j3;
	FREAL r1, r2, r3, r, drn, rd[3];
	thrust::complex<FREAL> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                    zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	thrust::complex<FREAL> zhi, zgi;

	float zhi_real, zhi_imag, zgi_real, zgi_imag;

	int iii, jjj;


#define zhelem(i, j, k) s_[3*num_warps*(i) + num_warps*(j) + (k)]
#define zgelem(i, j, k) (s_ + gelem_pad)[3*num_warps*(i) + num_warps*(j) + (k)]
	
	
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
	p[0][0] = -0.25*sm;
	p[1][0] =  0.25*sm;
	p[2][0] =  0.25*sp;
	p[3][0] = -0.25*sp;

	g1 = gi[ig];
	p1 = ome[ig];
	rp = 1 + g1;
	rm = 1 - g1;
	f[0] = 0.25*rm*sm;
	f[1] = 0.25*rp*sm;
	f[2] = 0.25*rp*sp;
	f[3] = 0.25*rm*sp;
	p[0][1] = -0.25*rm;
	p[1][1] = -0.25*rp;
	p[2][1] = 0.25*rp;
	p[3][1] = 0.25*rm;

	
   
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
	

	zwi = thrust::complex<FREAL>(0, fr);
	
	zc0 = ((FREAL) 1.0)/(((FREAL) 4.)*(pi)*(zge));
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

	zfhi    = (((FREAL) 1.) + ((FREAL) 1.)/zs2 - ((FREAL) 1.)/zzs)*zezs/r - zc2*(((FREAL) 1.)/zp2 - ((FREAL) 1.)/zzp)*zezp/r;
	zcappa  = (((FREAL) 1.) + ((FREAL) 3.)/zs2 - ((FREAL) 3.f)/zzs)*zezs/r - zc2*(((FREAL) 1.) + ((FREAL) 3.)/zp2 - ((FREAL) 3.)/zzp)*zezp/r;
	zfhidr  = (((FREAL) -2.)+ zzs + ((FREAL) 3.)/zzs - ((FREAL) 3.)/zs2)*zezs/(r*r) - zc2*(((FREAL) -1.) + ((FREAL) 3.)/zzp - ((FREAL) 3.)/zp2)*zezp/(r*r);
	zcappadr= (zzs - ((FREAL) 4.) + ((FREAL) 9.f)/zzs - ((FREAL) 9.)/zs2)*zezs/(r*r) - zc2*(zzp - ((FREAL) 4.) + ((FREAL) 9.)/zzp - ((FREAL) 9.f)/zp2)*zezp/(r*r);

	zaa = zfhidr-zcappa/r;
	zbb = ((FREAL) 4.)*zcappa/r -((FREAL) 2.)*zcappadr;
	zcc = (zc1-((FREAL) 2.))*(zaa + ((FREAL) 0.5)*zbb-((FREAL) 3.0)*zcappa/r)-((FREAL) 2.0)*zcappa/r;

	p12 = p1*p2*det;
	
    for (j = 0; j < 3; ++j)
    {	for (i = 0; i < 3; ++i)
        {
            zgi = (zc0*(zfhi*delta[j][i] - zcappa*rd[j]*rd[i]));
            

            zhi = (((FREAL) 1.0)/(((FREAL) 4.0)*pi))*((zaa*(drn*delta[j][i] + 
                                rd[j]*rn_cached[i])) + rd[i]*rd[j]*drn*zbb + 
                        rd[i]*rn_cached[j]*zcc);
        
            if (ii == jj && fast_singular)
            {
                zgi = zgi - (c1/r)*(c2*delta[j][i] + rd[i]*rd[j]);
                zhi = zhi - (c3/(r*r))*(drn*(c4*delta[j][i] + ((FREAL) 3.0)*rd[i]*rd[j]) + c4*(rd[j]*rn_cached[i] - rd[i]*rn_cached[j]));
            }
          
            zgi = zgi*p12;
			zgi_real = zgi.real();
			zgi_imag = zgi.imag();

			for (int offset = 16; offset > 0; offset = offset/2)
				zgi_real += __shfl_down(zgi_real, offset);
			for (int offset = 16; offset > 0; offset = offset/2)
				zgi_imag += __shfl_down(zgi_imag, offset);
			
            zhi = zhi*p12;
			zhi_real = zhi.real();
			zhi_imag = zhi.imag();
			for (int offset = 16; offset > 0; offset = offset/2)
				zhi_real += __shfl_down(zhi_real, offset);
			for (int offset = 16; offset > 0; offset = offset/2)
				zhi_imag += __shfl_down(zhi_imag, offset);

			if (lane == 0)
			{
				zhelem(j, i, warp) = thrust::complex<float>(zhi_real, zhi_imag);
				zgelem(j, i, warp) = thrust::complex<float>(zgi_real, zgi_imag);
			}
        }
    }
	__syncthreads();
	
	if (jg < 3 && ig < 3)
	{
		int index = 3*blockIdx.y + (3*nbe)*3*blockIdx.x + ig + (3*nbe)*jg;
		zg[index] = thrust::reduce(thrust::seq, &zgelem(jg, ig, 0), &zgelem(jg, ig, num_warps));
	} else if ((npg-3) <= jg && jg < npg && (npg-3) <= ig && ig < npg) //Split the warps
	{
		int index = 3*blockIdx.y + (3*nbe)*3*blockIdx.x + (ig-(npg-3)) + (3*nbe)*(jg-(npg-3));
		zh[index] = thrust::reduce(thrust::seq, &zhelem((jg-(npg-3)), (ig-(npg-3)), 0), &zhelem((jg-(npg-3)), (ig-(npg-3)), num_warps));
	}

//	} else if (3 <= jg && jg  < 6 && 3 <= ig && ig < 6)
//	{
//		int index = 3*blockIdx.y + (3*nbe)*3*blockIdx.x + (ig-3) + (3*nbe)*(jg-3);
//		zh[index] = thrust::reduce(thrust::seq, &zhelem((jg-3), (ig-3), 0), &zhelem((jg-3), (ig-3), npg*npg));
//	}
}


void cuda_ghmatecd_(int* nbe,
                    int* npg,
                    int* n,
                    int* np,
                    thrust::complex<FREAL>* zge,
                    thrust::complex<FREAL>* zcs,
                    thrust::complex<FREAL>* zcp,
                    FREAL* c1,
                    FREAL* c2,
                    FREAL* c3,
                    FREAL* c4,
                    FREAL* fr,
                    FREAL* zhest_,
                    FREAL* zgest_,
                    thrust::complex<FREAL>* zgp_,
                    thrust::complex<FREAL>* zhp_,
                    int* fast_singular,
                    int* status
                   )
{
	
	size_t column_size = 2*(3*(*nbe))*sizeof(thrust::complex<FREAL>);
	
	int shared_mem_size = 2*3*3*(*npg)*(*npg)*sizeof(thrust::complex<FREAL>);
	cudaError_t error;
	
	thrust::complex<FREAL>* device_zh;
	thrust::complex<FREAL>* device_zg;

	int* device_return_status;
	int return_status;

	/*Cast os parâmetros de volta para o tipo original*/
//	FREAL (*zhest)[3*(*nbe)] = (FREAL (*)[3*(*nbe)]) zhest_;
//	FREAL (*zgest)[3*(*nbe)] = (FREAL (*)[3*(*nbe)]) zgest_;
	thrust::complex<FREAL> (*zgp)[3*(*nbe)] = (thrust::complex<FREAL> (*)[3*(*nbe)]) zgp_;
	thrust::complex<FREAL> (*zhp)[3*(*nbe)] = (thrust::complex<FREAL> (*)[3*(*nbe)]) zhp_;

	int i, iterations, width;
	dim3 threadsPerBlock(*npg,*npg);

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	width = largest_possible_width(column_size, *nbe, &iterations);

	error = cudaMalloc(&device_zh, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

	error = cudaMalloc(&device_zg, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<FREAL>));
	cuda_assert(error);

	for (i = 0; i < iterations; ++i)
	{
		int starting_column = width*i;
//		if (starting_column + width > *n)
//			width = *n - starting_column;
		if (starting_column + width > *nbe)
			width = *nbe - starting_column;
		dim3 numBlocks(width, *nbe);


		error = cudaMemset(device_return_status, 0, sizeof(int));
		cuda_assert(error);

		error = cudaMemset(device_zh, 0, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<FREAL>));
		cuda_assert(error);

		error = cudaMemset(device_zg, 0, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<FREAL>));
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
							(FREAL (*)[3]) device_etas,
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
//							*npg,
							*n,
							*nbe,
							*n,
							starting_column,
                            *fast_singular,
							device_return_status
							);
		cudaDeviceSynchronize();

		error = cudaMemcpy(&return_status, device_return_status, sizeof(int), cudaMemcpyDeviceToHost);
		cuda_assert(error);

		if (return_status != 0)
		{
			fputs("Matriz Singular\n", stderr);
		}

		error = cudaMemcpy(&zhp[3*starting_column], device_zh, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<FREAL>), cudaMemcpyDeviceToHost);
		cuda_assert(error);
		error = cudaMemcpy(&zgp[3*starting_column], device_zg, (3*(*nbe))*(3*(width))*sizeof(thrust::complex<FREAL>), cudaMemcpyDeviceToHost);
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
