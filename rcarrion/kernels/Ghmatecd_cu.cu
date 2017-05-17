#include <cstdio>
#include <cmath>
#include <thrust/complex.h>

#define restrict __restrict__


extern "C"{

void cuda_assert(cudaError_t error)
{
	if (error != cudaSuccess)
	{
		fputs(cudaGetErrorString(cudaGetLastError()), stderr);
		exit(1);
	}
}


void sing_de_(thrust::complex<float> zhelem[3][3],
              thrust::complex<float> zgelem[3][3],
              float co[][4],
              float* cxm,
              float* cym,
              float* czm,
              float eta[],
              thrust::complex<float>* zge,
              thrust::complex<float>* zcs,
              thrust::complex<float>* zcp,
              float* c1,
              float* c2,
              float* c3,
              float* c4,
              float delta[][3],
              float* pi,
              float* fr,
			  float gi[],
			  float ome[],
              int* npg
              );

void nonsingd_(thrust::complex<float> zhelem[3][3],
               thrust::complex<float> zgelem[3][3],
               float co[][4],
               float* cxm,
               float* cym,
               float* czm,
               float eta[],
               thrust::complex<float>* zge,
               thrust::complex<float>* zcs,
               thrust::complex<float>* zcp,
               float delta[][3],
			   float* pi,
               float* fr,
               float gi[],
			   float ome[],
               int* npg
               );


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
						   float delta[3][3],
						   float pi,
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
						   int nx,
						   int ne,
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
    __shared__ float rn_cached[3];
	float g0, g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
    float cxg, cyg, czg, cxp, cyp, czp;
    float j1, j2, j3;
    float r1, r2, r3, r, drn, rd[3];
    thrust::complex<float> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                    zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	thrust::complex<float> zhi, zgi;

	__shared__ thrust::complex<float> zhelem[3][3], zgelem[3][3];

	//const float pi  = 3.141592654;
	
	int iii, jjj;

	if (threadIdx.x < 3 && threadIdx.y < 3)
	{
		zgelem[threadIdx.x][threadIdx.y] = 0;
		zhelem[threadIdx.x][threadIdx.x] = 0;
	}
	if (threadIdx.x < 4 && threadIdx.y == 1)
	{
		co[0][threadIdx.x] = cx[cone[ne*threadIdx.x + jj] - 1];
		co[1][threadIdx.x] = cy[cone[ne*threadIdx.x + jj] - 1];
		co[2][threadIdx.x] = cz[cone[ne*threadIdx.x + jj] - 1];
	}
	if (threadIdx.x < 3 && threadIdx.y == 1)
	{
		rn_cached[threadIdx.x] = rn[jj][threadIdx.x];
	}
	__syncthreads();

	cxp = cxm[ii];
	cyp = cym[ii];
	czp = czm[ii];

//	zh[index] = 0;
//	zg[index] = 0;
    
//	__syncthreads();

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
    
    /*
    if (threadIdx.y < 2 && threadIdx.x < 3)
    {
        int iii = threadIdx.y, jjj = threadIdx.x;
        xj[jjj][iii] = p[0][iii]*co[jjj][0] + p[1][iii]*co[jjj][1]+ p[2][iii]*co[jjj][2] + p[3][iii]*co[jjj][3];
    }
    __syncthreads();
    */

    j1 = xj[1][0]*xj[2][1]-xj[1][1]*xj[2][0];
	j2 = xj[0][1]*xj[2][0]-xj[0][0]*xj[2][1];
	j3 = xj[0][0]*xj[1][1]-xj[0][1]*xj[1][0];

	/*otimizar*/
	det = sqrt(j1*j1 + j2*j2 + j3*j3);
	/**/


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
	

	for (j = 0; j < 3; ++j)
	{	for (i = 0; i < 3; ++i)
		{
			int index = 3*ii + nx*3*jj + i + nx*j;
			zgi = (zc0*(zfhi*delta[j][i] - zcappa*rd[j]*rd[i]));
			

			zhi = (1.0f/(4.0f*pi))*((zaa*(drn*delta[j][i] + 
								rd[j]*rn_cached[i])) + rd[i]*rd[j]*drn*zbb + 
						rd[i]*rn_cached[j]*zcc);
		
            
			if (ii == jj)
			{
				zgi = zgi - (c1/r)*(c2*delta[j][i] + rd[i]*rd[j]);
				zhi = zhi - (c3/(r*r))*(drn*(c4*delta[j][i] + 3.0f*rd[i]*rd[j]) + c4*(rd[j]*rn_cached[i] - rd[i]*rn_cached[j]));
			}
			
			zgi = zgi*p12;
			zhi = zhi*p12;

			atomicAdd(((float*) &zgelem[j][i])    , zgi.real());
			atomicAdd(((float*) &zgelem[j][i]) + 1, zgi.imag());

			atomicAdd(((float*) &zhelem[j][i])    , zhi.real());
			atomicAdd(((float*) &zhelem[j][i]) + 1, zhi.imag());

//			atomicAdd(((float*) &zg[index])    , zgi.real());
//			atomicAdd(((float*) &zg[index]) + 1, zgi.imag());
	
//			atomicAdd(((float*) &zh[index])    , zhi.real());
//			atomicAdd(((float*) &zh[index]) + 1, zhi.imag());
		
		}
	}

	__syncthreads();
	
	if (jg < 3 && ig < 3)
	{
		int index = 3*ii + nx*3*jj + ig + nx*jg;
		zg[index] = zgelem[jg][ig];
		zh[index] = zhelem[jg][ig];
	}
}


__global__ void nonsingd_kernel(
						   thrust::complex<float> zhelem[3][3],
						   thrust::complex<float> zgelem[3][3],
						   float co[3][4],
						   float cxp,
						   float cyp,
						   float czp,
						   float rn[3],
						   thrust::complex<float> zge,
						   thrust::complex<float> zcs,
						   thrust::complex<float> zcp,
						   float delta[3][3],
						   float pi,
						   float fr,
						   float gi[],
						   float ome[],
						   int npg
						   )
{

	const int ig = blockIdx.y;
	const int jg = blockIdx.x;


	float p[4][2], xj[3][2], f[4];
	float g0, g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
    float cxg, cyg, czg; //, cxp, cyp, czp;
    float j1, j2, j3;
    float r1, r2, r3, r, drn, rd[3];
    thrust::complex<float> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                    zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	thrust::complex<float> zhi, zgi;

	//const float pi  = 3.141592654;
	
	int ii, jj;
/*
	n1 = cone[ne*0 + j];
	n2 = cone[ne*1 + j];
	n3 = cone[ne*2 + j];
	n4 = cone[ne*3 + j];

	co[0][0] = cx[n1 - 1];
	co[1][0] = cy[n1 - 1];
	co[2][0] = cz[n1 - 1];
	co[0][1] = cx[n2 - 1];
	co[1][1] = cy[n2 - 1];
	co[2][1] = cz[n2 - 1];
	co[0][2] = cx[n3 - 1];
	co[1][2] = cy[n3 - 1];
	co[2][2] = cz[n3 - 1];
	co[0][3] = cx[n4 - 1];
	co[1][3] = cy[n4 - 1];
	co[2][3] = cz[n4 - 1];

	cxp = cxm[i];
	cyp = cym[i];
	czp = czm[i];
*/
	zhelem[threadIdx.x][threadIdx.y] = 0;
	zgelem[threadIdx.x][threadIdx.y] = 0;
    
	__syncthreads();

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

	for (ii = 0; ii < 2; ++ii)
	{
		for (jj = 0; jj < 3; ++jj)
		{
			xj[jj][ii] = p[0][ii]*co[jj][0] + p[1][ii]*co[jj][1]+ p[2][ii]*co[jj][2] + p[3][ii]*co[jj][3];
		}
	}

    j1 = xj[1][0]*xj[2][1]-xj[1][1]*xj[2][0];
	j2 = xj[0][1]*xj[2][0]-xj[0][0]*xj[2][1];
	j3 = xj[0][0]*xj[1][1]-xj[0][1]*xj[1][0];

	/*otimizar*/
	det = sqrt(j1*j1 + j2*j2 + j3*j3);
	/**/

/*
	if (det < 1e-5)
	{
		//algo deveria acontecer
	}
*/

	cxg = 0;
	cyg = 0;
	czg = 0;
	for (ii = 0; ii < 4; ++ii)
	{
		cxg = cxg + co[0][ii]*f[ii];
		cyg = cyg + co[1][ii]*f[ii];
		czg = czg + co[2][ii]*f[ii];
	}


	r1    = cxg - cxp;
	r2    = cyg - cyp;
	r3    = czg - czp;

	/*otimizar*/
	r     = sqrt(r1*r1 + r2*r2 + r3*r3);
	drn   = (r1*rn[0] + r2*rn[1] + r3*rn[2])/r;
	rd[0] = r1/r;
	rd[1] = r2/r;
	rd[2] = r3/r;
	/**/

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

	zgi = (zc0*(zfhi*delta[threadIdx.x][threadIdx.y] - zcappa*rd[threadIdx.x]*rd[threadIdx.y]));
	zhi = (1.0f/(4.0f*pi))*((zaa*(drn*delta[threadIdx.x][threadIdx.y] + 
						rd[threadIdx.x]*rn[threadIdx.y])) + rd[threadIdx.y]*rd[threadIdx.x]*drn*zbb + 
				rd[threadIdx.y]*rn[threadIdx.x]*zcc);
	zgi = zgi*p12;
	zhi = zhi*p12;

	atomicAdd((float*) &zgelem[threadIdx.x][threadIdx.y]    , zgi.real());
	atomicAdd((float*) &zgelem[threadIdx.x][threadIdx.y] + 1, zgi.imag());

	atomicAdd((float*) &zhelem[threadIdx.x][threadIdx.y]    , zhi.real());
	atomicAdd((float*) &zhelem[threadIdx.x][threadIdx.y] + 1, zhi.imag());
}
void cuda_ghmatecd_(int* ne,
                    int* nbe,
                    int* nx,
                    int* npg,
                    int* ncox,
					int* n,
					int* cone_,
					float cx[],
                    float cy[],
                    float cz[],
                    float cxm[],
                    float cym[],
                    float czm[],
                    float etas[][3],
                    thrust::complex<float>* zge,
                    thrust::complex<float>* zcs,
                    thrust::complex<float>* zcp,
                    float* c1,
                    float* c2,
                    float* c3,
                    float* c4,
                    float delta[3][3],
                    float* fr,
                    float* zhest_,
					float* zgest_,
					thrust::complex<float>* zgp_,
                    thrust::complex<float>* zhp_,
					float ome[],
					float* gi,
					int* status
                   )
{
	dim3 threadsPerBlock(*npg,*npg);
	dim3 numBlocks(*n, *nbe);
	cudaError_t error;
    
    int n1, n2, n3, n4;
    int i, j, ii, jj;
	thrust::complex<float> zhelem[3][3];
	thrust::complex<float> zgelem[3][3];
	float co[3][4];

	thrust::complex<float>* device_zh;
	thrust::complex<float>* device_zg;
	float* device_cx;
	float* device_cy;
	float* device_cz;
	float* device_cxm;
	float* device_cym;
	float* device_czm;
	float* device_gi;
	float* device_ome;
	float* device_etas;
	float* device_delta;
	int* device_cone;

	int* device_return_status;
	int return_status;

	/*Cast os parâmetros de volta para o tipo original*/
	int (*cone)[*ne]          = (int (*)[*ne])   cone_;
	float (*zhest)[*nx]       = (float (*)[*nx]) zhest_;
	float (*zgest)[*nx]       = (float (*)[*nx]) zgest_;
	thrust::complex<float> (*zgp)[*nx] = (thrust::complex<float> (*)[*nx]) zgp_;
	thrust::complex<float> (*zhp)[*nx] = (thrust::complex<float> (*)[*nx]) zhp_;

	float pi  = 3.141592654;

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	error = cudaMalloc(&device_cone, 4*(*ne)*sizeof(int));
	cuda_assert(error);

	error = cudaMalloc(&device_zh, (*nx)*(*nx)*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMalloc(&device_zg, (*nx)*(*nx)*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMemset(device_return_status, 0, sizeof(int));
	cuda_assert(error);

	error = cudaMemset(device_zh, 0, (*nx)*(*nx)*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMemset(device_zg, 0, (*nx)*(*nx)*sizeof(thrust::complex<float>));
	cuda_assert(error);

	error = cudaMalloc(&device_cx, (*ncox)*sizeof(float));
	cuda_assert(error);

	error = cudaMalloc(&device_cy, (*ncox)*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_cz, (*ncox)*sizeof(float));
	cuda_assert(error);

	error = cudaMalloc(&device_cxm, (*ne)*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_cym, (*ne)*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_czm, (*ne)*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_gi, (*npg)*sizeof(float));
	cuda_assert(error);

	error = cudaMalloc(&device_ome, (*npg)*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_etas, (*nx)*3*sizeof(float));
	cuda_assert(error);
	
	error = cudaMalloc(&device_delta, 3*3*sizeof(float));
	cuda_assert(error);

	error = cudaMemcpy(device_cone, cone, 4*(*ne)*sizeof(int), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_cx, cx, (*ncox)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_cy, cy, (*ncox)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_cz, cz, (*ncox)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_cxm, cxm, (*ne)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_cym, cym, (*ne)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_czm, czm, (*ne)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_gi, gi, (*npg)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_ome, ome, (*npg)*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_etas, etas, (*n)*3*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_delta, delta, 3*3*sizeof(float), cudaMemcpyHostToDevice);
	cuda_assert(error);

	cudaDeviceSynchronize();

	printf("n = %d, nbe = %d, npg = %d\n", *n,*nbe, *npg);

//	cudaFuncSetCacheConfig(ghmatecd_kernel, cudaFuncCachePreferShared);
	ghmatecd_kernel<<<numBlocks, threadsPerBlock>>>(
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
						(float (*)[3]) device_delta,
						pi,
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
						*nx,
						*ne,
						device_return_status
						);
	cudaDeviceSynchronize();

	error = cudaMemcpy(&return_status, device_return_status, sizeof(int), cudaMemcpyDeviceToHost);
	cuda_assert(error);

	if (return_status != 0)
	{
		fputs("Matriz Singular", stderr);
	}

	error = cudaMemcpy(zhp_, device_zh, (*nx)*(*nx)*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
	cuda_assert(error);
	error = cudaMemcpy(zgp_, device_zg, (*nx)*(*nx)*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
	cuda_assert(error);

	
	
	
	for (i = 0; i < *nbe; ++i)
	{
		ii = 3*i;
		for (int jjj = 0; jjj < 3; ++jjj)
		{   for (int iii = 0; iii < 3; ++iii)
			{
				zgp[jjj+ii][iii+ii] += zgest[jjj+ii][iii+ii];
				zhp[jjj+ii][iii+ii] += zhest[jjj+ii][iii+ii];
			}
		}
	}

/*
	for (j = 0; j < *n; ++j)
    {
		n1 = cone[0][j];
        n2 = cone[1][j];
        n3 = cone[2][j];
        n4 = cone[3][j];

        co[0][0] = cx[n1 - 1];
        co[1][0] = cy[n1 - 1];
        co[2][0] = cz[n1 - 1];
        co[0][1] = cx[n2 - 1];
        co[1][1] = cy[n2 - 1];
        co[2][1] = cz[n2 - 1];
        co[0][2] = cx[n3 - 1];
        co[1][2] = cy[n3 - 1];
        co[2][2] = cz[n3 - 1];
        co[0][3] = cx[n4 - 1];
        co[1][3] = cy[n4 - 1];
        co[2][3] = cz[n4 - 1];

//		error = cudaMemcpy(device_co, co,  3*4*sizeof(float), cudaMemcpyHostToDevice);
//		cuda_assert(error);
       		
		jj = 3*j;
        for (i = 0; i < *nbe; ++i)
        {
            ii = 3*i;
	//		printf("i = %d, j = %d\n", i, j);
            if (i == j)
            {

				sing_de_(zhelem, 
                         zgelem, 
                         co, 
                         &cxm[i], 
                         &cym[i], 
                         &czm[i], 
                         etas[j],
                         zge,
                         zcs,
                         zcp,
                         c1,
                         c2,
                         c3,
                         c4,
                         delta,
                         &pi,
                         fr,
						 gi,
						 ome,
                         npg
                        );
                
				for (int jjj = 0; jjj < 3; ++jjj)
				{   for (int iii = 0; iii < 3; ++iii)
                    {
                    
						printf("error: %f\n", abs(zgp[jjj+jj][iii+ii] - zgelem[jjj][iii]));
						if (abs(zgp[jjj+jj][iii+ii] - zgelem[jjj][iii]) > 1E-5)
						{
							printf("error: %f\n", abs(zgp[jjj+jj][iii+ii] - zgelem[jjj][iii]));
						}
						
						zgp[jjj+jj][iii+ii] = zgelem[jjj][iii] + zgest[jjj+jj][iii+ii];
                        zhp[jjj+jj][iii+ii] = zhelem[jjj][iii] + zhest[jjj+jj][iii+ii];
					}
				}
            }

//			else
//          {
//		
//                nonsingd_kernel<<<numBlocks, threadsPerBlock>>>(
//						(thrust::complex<float> (*)[3]) device_zhelem,
//						(thrust::complex<float> (*)[3]) device_zgelem,
//						(float (*)[4])          device_co,
//						cxm[i],
//						cym[i],
//						czm[i],
//						&device_etas[3*j],
//						*zge,
//						*zcs,
//						*zcp,
//						(float (*)[3])device_delta,
//						pi,
//						*fr,
//						device_gi,
//						device_ome,
//						*npg
//						);
//				cudaDeviceSynchronize();				

                
//                nonsingd_(zhelem, 
//                          zgelem, 
//                          co, 
//                          &cxm[i], 
//                          &cym[i], 
//                          &czm[i],
//                          etas[j],
//                          zge,
//                          zcs,
//                          zcp,
//                          delta,
//                          &pi,
//                          fr,
//						  gi,
//						  ome,
//                          npg
//                         );
//					
				
//				error = cudaMemcpy(zgelem, device_zgelem, 3*3*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
//				cuda_assert(error);
//				error = cudaMemcpy(zhelem, device_zhelem, 3*3*sizeof(thrust::complex<float>), cudaMemcpyDeviceToHost);
//				cuda_assert(error);

				

//				for (int jjj = 0; jjj < 3; ++jjj)
//				{   for (int iii = 0; iii < 3; ++iii)
//                    {
//						if (abs(zgp[jjj+jj][iii+ii] - zgelem[jjj][iii]) > 1E-3)
//						{
//							printf("zgp[jjj+jj][iii+ii] = %f,%f\n", zgp[jjj+jj][iii+ii].real(), zgp[jjj+jj][iii+ii].imag());
//						}
//						zgp[jjj+jj][iii+ii] = zgelem[jjj][iii];
//                        zhp[jjj+jj][iii+ii] = zhelem[jjj][iii];
//                    }
//				}
//			}
			
        }
    }
	
*/	
	error = cudaFree(device_cone);
	cuda_assert(error);
	error = cudaFree(device_zh);
	cuda_assert(error);
	error = cudaFree(device_zg);
	cuda_assert(error);
	error = cudaFree(device_gi);
	cuda_assert(error);
	error = cudaFree(device_ome);
	cuda_assert(error);
	error = cudaFree(device_etas);
	cuda_assert(error);
	error = cudaFree(device_delta);
	cuda_assert(error);
//	error = cudaFree(device_co);
//	cuda_assert(error);
	error = cudaFree(device_cx);
	cuda_assert(error);
	error = cudaFree(device_cy);
	cuda_assert(error);
	error = cudaFree(device_cz);
	cuda_assert(error);
	error = cudaFree(device_cxm);
	cuda_assert(error);
	error = cudaFree(device_cym);
	cuda_assert(error);
	error = cudaFree(device_czm);
	cuda_assert(error);
}
}
