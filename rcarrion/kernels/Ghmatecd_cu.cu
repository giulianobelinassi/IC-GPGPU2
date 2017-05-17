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

void build_cos(
		int npg, 
		int ne, 
		float output_cos[][3][4], 
		int* cone_, 
		float cx[], 
		float cy[], 
		float cz[]
		)
{
	int j, i;
	int (*cone)[ne] = (int (*)[ne]) cone_;
	int n_[4];	

	for (j = 0; j < npg; ++j)
	{
		n_[0] = cone_[ne*0 + j];
		n_[1] = cone_[ne*1 + j];
		n_[2] = cone_[ne*2 + j];
		n_[3] = cone_[ne*3 + j];

		output_cos[j][0][0] = cx[n_[0] - 1];
		output_cos[j][1][0] = cy[n_[0] - 1];
		output_cos[j][2][0] = cz[n_[0] - 1];
		output_cos[j][0][1] = cx[n_[1] - 1];
		output_cos[j][1][1] = cy[n_[1] - 1];
		output_cos[j][2][1] = cz[n_[1] - 1];
		output_cos[j][0][2] = cx[n_[2] - 1];
		output_cos[j][1][2] = cy[n_[2] - 1];
		output_cos[j][2][2] = cz[n_[2] - 1];
		output_cos[j][0][3] = cx[n_[3] - 1];
		output_cos[j][1][3] = cy[n_[3] - 1];
		output_cos[j][2][3] = cz[n_[3] - 1];
	}
}

__device__ thrust::complex<float> reduct(thrust::complex<float> data[], int length)
{	
	int i;
	thrust::complex<float> sum = 0;
	for (i = 0; i < length; ++i)
		sum += data[i];
	return sum;
}

__global__ void build_3x3_matrix(
		thrust::complex<float> zh[],
		thrust::complex<float> zg[],
		thrust::complex<float> zc0,
		thrust::complex<float> zfhi,
		float delta[3][3],
		thrust::complex<float> zcappa,
		float r1,
		float r2,
		float r3,
		float rn[3],
		float drn,
		thrust::complex<float> zaa,
		thrust::complex<float> zbb,
		thrust::complex<float> zcc,
		float r,
		float c1,
		float c2,
		float c3,
		float c4,
		float pi,
		float p12,
		int nx,
		int ii,
		int jj
		)
{
	const int i = threadIdx.y;
	const int j = threadIdx.x;
	float rd[3] = {r1/r, r2/r, r3/r};
	int index = 3*ii + nx*3*jj + i + nx*j;
	
	thrust::complex<float> zgi, zhi;
	zgi = (zc0*(zfhi*delta[j][i] - zcappa*rd[j]*rd[i]));
	
	zhi = (1.0f/(4.0f*pi))*((zaa*(drn*delta[j][i] + 
						rd[j]*rn[i])) + rd[i]*rd[j]*drn*zbb + 
				rd[i]*rn[j]*zcc);

	if (ii == jj)
	{
		zgi = zgi - (c1/r)*(c2*delta[j][i] + rd[i]*rd[j]);
		zhi = zhi - (c3/(r*r))*(drn*(c4*delta[j][i] + 3.0f*rd[i]*rd[j]) + c4*(rd[j]*rn[i] - rd[i]*rn[j]));
	}
	
	zgi = zgi*p12;
	zhi = zhi*p12;

	atomicAdd(((float*) &zg[index])    , zgi.real());
	atomicAdd(((float*) &zg[index]) + 1, zgi.imag());

	atomicAdd(((float*) &zh[index])    , zhi.real());
	atomicAdd(((float*) &zh[index]) + 1, zhi.imag());
}


__global__ void ghmatecd_kernel(
	float co[][3][4],
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

	dim3 threads3x3(3,3);

	float p[4][2], f[4];
	float xj[3][2];
	float g0, g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
	float cxg, cyg, czg, cxp, cyp, czp;
	float j1, j2, j3;
	float r1, r2, r3, r, drn, rd[3];
	thrust::complex<float> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs,
                               zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	thrust::complex<float> zhi, zgi;

//	__shared__ thrust::complex<float> zgelem[3][3][64];
//	__shared__ thrust::complex<float> zhelem[3][3][64];
	//const float pi  = 3.141592654;
	
	int iii, jjj;

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
			xj[jjj][iii] = p[0][iii]*co[jj][jjj][0] + p[1][iii]*co[jj][jjj][1]+ p[2][iii]*co[jj][jjj][2] + p[3][iii]*co[jj][jjj][3];
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
		cxg = cxg + co[jj][0][iii]*f[iii];
		cyg = cyg + co[jj][1][iii]*f[iii];
		czg = czg + co[jj][2][iii]*f[iii];
	}

	r1    = cxg - cxp;
	r2    = cyg - cyp;
	r3    = czg - czp;
	
	r     = sqrt(r1*r1 + r2*r2 + r3*r3);
	drn   = (r1*rn[jj][0] + r2*rn[jj][1] + r3*rn[jj][2])/r;
/*
	rd[0] = r1/r;
	rd[1] = r2/r;
	rd[2] = r3/r;
*/	
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
	
	build_3x3_matrix<<<threads3x3, 1>>>(
			zh,
			zg,
			zc0,
			zfhi,
			delta,
			zcappa,
			r1,
			r2,
			r3,
			rn[jj],
			drn,
			zaa,
			zbb,
			zcc,
			r,
			c1,
			c2,
			c3,
			c4,
			pi,
			p12,
			nx,
			ii,
			jj
	);
	/*
	for (j = 0; j < 3; ++j)
	{	for (i = 0; i < 3; ++i)
		{
			int index = 3*ii + nx*3*jj + i + nx*j;
			zgi = (zc0*(zfhi*delta[j][i] - zcappa*rd[j]*rd[i]));
			

			zhi = (1.0f/(4.0f*pi))*((zaa*(drn*delta[j][i] + 
								rd[j]*rn[jj][i])) + rd[i]*rd[j]*drn*zbb + 
						rd[i]*rn[jj][j]*zcc);
		
            
			if (ii == jj)
			{
				zgi = zgi - (c1/r)*(c2*delta[j][i] + rd[i]*rd[j]);
				zhi = zhi - (c3/(r*r))*(drn*(c4*delta[j][i] + 3.0f*rd[i]*rd[j]) + c4*(rd[j]*rn[jj][i] - rd[i]*rn[jj][j]));
			}
			
			zgi = zgi*p12;
			zhi = zhi*p12;

//			zgelem[j][i][npg*jg + ig] = zgi;
//			zhelem[j][i][npg*jg + ig] = zhi;


			atomicAdd(((float*) &zg[index])    , zgi.real());
			atomicAdd(((float*) &zg[index]) + 1, zgi.imag());
	
			atomicAdd(((float*) &zh[index])    , zhi.real());
			atomicAdd(((float*) &zh[index]) + 1, zhi.imag());
		
		}
	}
	*/
/*
	__syncthreads();
	
	if (threadIdx.x < 3 && threadIdx.y < 3)
	{	
		int index = 3*ii + nx*3*jj + ig + nx*jg;
		zg[index] = reduct(zgelem[ig][jg], 64);
		zh[index] = reduct(zhelem[ig][jg], 64);
	}
	__syncthreads();
*/	
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
//	thrust::complex<float> zhelem[3][3];
//	thrust::complex<float> zgelem[3][3];
	float co[*npg][3][4];

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
	float* device_cos;
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

	build_cos(*npg, *ne, co, cone_, cx, cy, cz);

	error = cudaMalloc(&device_return_status, sizeof(int));
	cuda_assert(error);

	error = cudaMalloc(&device_cos, *npg*3*4*sizeof(float));
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

	error = cudaMemcpy(device_cos, co, 3*4*(*npg)*sizeof(float), cudaMemcpyHostToDevice);
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

	ghmatecd_kernel<<<numBlocks, threadsPerBlock>>>(
						(float (*)[3][4])device_cos,
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
	error = cudaFree(device_cos);
	cuda_assert(error);

}
}
