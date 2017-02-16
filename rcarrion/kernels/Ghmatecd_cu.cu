#include <cstdio>
#include <cmath>
#include "cuda_complex.hpp"

#define restrict __restrict__


extern "C"{

__device__ double myAtomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
		 unsigned long long int old = *address_as_ull, assumed;
		 do {
		     assumed = old;
			 old = atomicCAS(address_as_ull, assumed,
			 __double_as_longlong(val + __longlong_as_double(assumed)));
	     } while (assumed != old);
	    return __longlong_as_double(old);
}

void cuda_assert(cudaError_t error)
{
	if (error != cudaSuccess)
	{
		fputs(cudaGetErrorString(cudaGetLastError()), stderr);
		exit(1);
	}
}


void sing_de_(complex<double> zhelem[3][3],
              complex<double> zgelem[3][3],
              double co[][4],
              double* cxm,
              double* cym,
              double* czm,
              double eta[],
              complex<double>* zge,
              complex<double>* zcs,
              complex<double>* zcp,
              double* c1,
              double* c2,
              double* c3,
              double* c4,
              double delta[][3],
              double* pi,
              double* fr,
			  double gi[],
			  double ome[],
              int* npg
              );

void nonsingd_(complex<double> zhelem[3][3],
               complex<double> zgelem[3][3],
               double co[][4],
               double* cxm,
               double* cym,
               double* czm,
               double eta[],
               complex<double>* zge,
               complex<double>* zcs,
               complex<double>* zcp,
               double delta[][3],
               double* pi,
               double* fr,
               int* npg
               );


__global__ void nonsingd_kernel(
						   complex<float> zhelem[3][3],
						   complex<float> zgelem[3][3],
						   double co[3][4],
						   double cxp,
						   double cyp,
						   double czp,
						   double rn[3],
						   complex<double> zge,
						   complex<double> zcs,
						   complex<double> zcp,
						   double delta[3][3],
						   double pi,
						   double fr,
						   double gi[],
						   double ome[],
						   int npg
						   )
{

	const int ig = blockIdx.y;
	const int jg = blockIdx.x;


	double p[4][2], xj[3][2], f[4];
	double g0, g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
    double cxg, cyg, czg; //, cxp, cyp, czp;
    double j1, j2, j3;
    double r1, r2, r3, r, drn, rd[3];
    complex<double> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                    zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	complex<float> zhi, zgi;

	//const double pi  = 3.141592654;
	
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

	det = sqrt(j1*j1 + j2*j2 + j3*j3);

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
	r     = sqrt(r1*r1 + r2*r2 + r3*r3);
	drn   = (r1*rn[0] + r2*rn[1] + r3*rn[2])/r;
	rd[0] = r1/r;
	rd[1] = r2/r;
	rd[2] = r3/r;

	zwi = complex<double>(0, fr);
	
	zc0 = 1.0/(4*(pi)*(zge));
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

	zfhi    = (1. + 1./zs2 - 1./zzs)*zezs/r - zc2*(1./zp2 - 1./zzp)*zezp/r;
	zcappa  = (1. + 3./zs2 - 3./zzs)*zezs/r - zc2*(1. + 3./zp2 - 3./zzp)*zezp/r;
	zfhidr  = (-2.+ zzs + 3./zzs - 3./zs2)*zezs/(r*r) - zc2*(-1. + 3./zzp - 3./zp2)*zezp/(r*r);
	zcappadr= (zzs - 4. + 9./zzs - 9./zs2)*zezs/(r*r) - zc2*(zzp - 4. + 9./zzp - 9./zp2)*zezp/(r*r);

	zaa = zfhidr-zcappa/r;
	zbb = 4.*zcappa/r - 2.*zcappadr;
	zcc = (zc1-2.)*(zaa + 0.5*zbb-3.0*zcappa/r)-2.0*zcappa/r;

	p12 = p1*p2*det;

	zgi = (zc0*(zfhi*delta[threadIdx.x][threadIdx.y] - zcappa*rd[threadIdx.x]*rd[threadIdx.y])*zzp)*p12;
	zhi = (1.0/(4.0*pi))*((zaa*(drn*delta[threadIdx.x][threadIdx.y] + 
						rd[threadIdx.x]*rn[threadIdx.y])) + rd[threadIdx.y]*rd[threadIdx.x]*drn*zbb + 
				rd[threadIdx.y]*rn[threadIdx.x]*zcc)*p12;

	atomicAdd((float*) &zgelem[threadIdx.x][threadIdx.y]                , real(zgi));
	atomicAdd((float*) &zgelem[threadIdx.x][threadIdx.y] + sizeof(float), imag(zgi));

	atomicAdd((float*) &zhelem[threadIdx.x][threadIdx.x]                , real(zhi));
	atomicAdd((float*) &zhelem[threadIdx.x][threadIdx.x] + sizeof(float), imag(zhi));
}


__global__ void ghmatecd_kernel(int ne,
		                   int nbe,
						   int nx,
						   int npg,
						   int ncox,
						   int n,
						   complex<double>* zh,
						   complex<double>* zg,
						   complex<double> zge,
						   complex<double> zcs,
						   complex<double> zcp,
						   double* zhest_,
						   double* zgest_,
						   double cx[],
						   double cy[],
						   double cz[],
						   double* cxm,
						   double* cym,
						   double* czm,
						   double gi[],
						   double ome[],
						   double delta[3][3],
						   double rn[][3],
						   double fr,
						   int* cone
						   )
{
	const int i = blockIdx.y;
	const int j = blockIdx.x;
	const int index = 3*blockIdx.y + 3*blockIdx.x*nx + threadIdx.x*nx + threadIdx.y;
	
	__shared__ double p[4][2], xj[3][2], f[4];
	double g0, g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
    double cxg, cyg, czg, cxp, cyp, czp;
    double j1, j2, j3;
	double co[3][4];
    double r1, r2, r3, r, drn, rd[3];
    complex<double> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                    zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	
	const double pi  = 3.141592654;
	
	__shared__ int n1;
	__shared__ int n2;
	__shared__ int n3;
	__shared__ int n4;

	n1 = cone[ne*0 + j];
	n2 = cone[ne*1 + j];
	n3 = cone[ne*2 + j];
	n4 = cone[ne*3 + j];

	int ig, jg, ii, jj;

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

	zh[index] = 0;
	zg[index] = 0;
    
	for (jg = 0; jg < npg; ++jg)
	{
	
		g2 = gi[jg];
		p2 = ome[jg];
		sp = 1 + g2;
		sm = 1 - g2;
		p[0][0] = -0.25*sm;
		p[1][0] =  0.25*sm;
		p[2][0] =  0.25*sp;
		p[3][0] = -0.25*sp;

		for (ig = 0; ig < npg; ++ig)
		{
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

			__syncthreads();

			xj[threadIdx.x][threadIdx.y] = p[0][threadIdx.y]*co[threadIdx.x][0] + p[1][threadIdx.y]*co[threadIdx.x][1]+ p[2][threadIdx.y]*co[threadIdx.x][2] + p[3][threadIdx.y]*co[threadIdx.x][3];

			__syncthreads();
			/*
            for (ii = 0; ii < 2; ++ii)
            {
                for (jj = 0; jj < 3; ++jj)
                {
                    xj[jj][ii] = p[0][ii]*co[jj][0] + p[1][ii]*co[jj][1]+ p[2][ii]*co[jj][2] + p[3][ii]*co[jj][3];
                }
            }
*/
            j1 = xj[1][0]*xj[2][1]-xj[1][1]*xj[2][0];
            j2 = xj[0][1]*xj[2][0]-xj[0][0]*xj[2][1];
            j3 = xj[0][0]*xj[1][1]-xj[0][1]*xj[1][0];

            det = sqrt(j1*j1 + j2*j2 + j3*j3);

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
			r     = sqrt(r1*r1 + r2*r2 + r3*r3);
			drn   = (r1*rn[j][0] + r2*rn[j][1] + r3*rn[j][2])/r;
			rd[0] = r1/r;
			rd[1] = r2/r;
			rd[2] = r3/r;

			zwi = complex<double>(0, fr);
			
			zc0 = 1.0/(4*(pi)*(zge));
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

			zfhi    = (1. + 1./zs2 - 1./zzs)*zezs/r - zc2*(1./zp2 - 1./zzp)*zezp/r;
			zcappa  = (1. + 3./zs2 - 3./zzs)*zezs/r - zc2*(1. + 3./zp2 - 3./zzp)*zezp/r;
			zfhidr  = (-2.+ zzs + 3./zzs - 3./zs2)*zezs/(r*r) - zc2*(-1. + 3./zzp - 3./zp2)*zezp/(r*r);
			zcappadr= (zzs - 4. + 9./zzs - 9./zs2)*zezs/(r*r) - zc2*(zzp - 4. + 9./zzp - 9./zp2)*zezp/(r*r);

			zaa = zfhidr-zcappa/r;
			zbb = 4.*zcappa/r - 2.*zcappadr;
			zcc = (zc1-2.)*(zaa + 0.5*zbb-3.0*zcappa/r)-2.0*zcappa/r;

			p12 = p1*p2*det;

			zg[index] += (zc0*(zfhi*delta[threadIdx.x][threadIdx.y] - zcappa*rd[threadIdx.x]*rd[threadIdx.y]))*p12;
		    zh[index] += ((1.0/(4.0*pi))*((zaa*(drn*delta[threadIdx.x][threadIdx.y] + rd[threadIdx.x]*rn[blockIdx.x][threadIdx.y])) + rd[threadIdx.y]*rd[threadIdx.x]*drn*zbb + rd[threadIdx.y]*rn[blockIdx.x][threadIdx.x]*zcc))*p12;
		}
	}
}

void cuda_ghmatecd_(int* ne,
                    int* nbe,
                    int* nx,
                    int* npg,
                    int* ncox,
					int* n,
					int* cone_,
					double cx[],
                    double cy[],
                    double cz[],
                    double cxm[],
                    double cym[],
                    double czm[],
                    double etas[][3],
                    complex<double>* zge,
                    complex<double>* zcs,
                    complex<double>* zcp,
                    double* c1,
                    double* c2,
                    double* c3,
                    double* c4,
                    double delta[3][3],
                    double* fr,
                    double* zhest_,
					double* zgest_,
					complex<double>* zgp_,
                    complex<double>* zhp_,
					double ome[],
					double* gi,
					int* status
                   )
{
	dim3 threadsPerBlock(3,3);
	dim3 numBlocks(*n, *nbe);
	cudaError_t error;

    double co[3][4];
    int n1, n2, n3, n4;
    int i, j, ii, jj;
    complex<double> zgelem[3][3];
	complex<double> zhelem[3][3];
	complex<float> zgelemf[3][3];
	complex<float> zhelemf[3][3];

	complex<float>* device_zhelem;
	complex<float>* device_zgelem;
	double* device_gi;
	double* device_ome;
	double* device_etas;
	double* device_delta;
	double* device_co;


	/*Cast os parâmetros de volta para o tipo original*/
	int (*cone)[*ne]           = (int (*)[*ne]) cone_;
	double (*zhest)[*nx]       = (double (*)[*nx]) zhest_;
	double (*zgest)[*nx]       = (double (*)[*nx]) zgest_;
	complex<double> (*zgp)[*nx] = (complex<double> (*)[*nx]) zgp_;
	complex<double> (*zhp)[*nx] = (complex<double> (*)[*nx]) zhp_;

	double pi  = 3.141592654;

	error = cudaMalloc(&device_zhelem, 3*3*sizeof(complex<float>));
	cuda_assert(error);

	error = cudaMalloc(&device_zgelem, 3*3*sizeof(complex<float>));
	cuda_assert(error);

	error = cudaMalloc(&device_gi, (*npg)*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_ome, (*npg)*sizeof(double));
	cuda_assert(error);
	
	error = cudaMalloc(&device_etas, (*nx)*3*sizeof(double));
	cuda_assert(error);
	
	error = cudaMalloc(&device_delta, 3*3*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_co,  3*4*sizeof(double));
	cuda_assert(error);
	
	error = cudaMemcpy(device_gi, gi, (*npg)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_ome, ome, (*npg)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_etas, etas, (*n)*3*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);
	

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

		error = cudaMemcpy(device_co, (double*) co,  3*4*sizeof(double), cudaMemcpyHostToDevice);
		cuda_assert(error);
        
		jj = 3*j;
        for (i = 0; i < *nbe; ++i)
        {
            ii = 3*i;

			printf("i = %d, j = %d\n", i, j);
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
                        zgp[jjj+jj][iii+ii] = zgelem[jjj][iii] + zgest[jjj+jj][iii+ii];
                        zhp[jjj+jj][iii+ii] = zhelem[jjj][iii] + zhest[jjj+jj][iii+ii];
                    }
				}
            }
            else
            {
				nonsingd_kernel<<<numBlocks, threadsPerBlock>>>(
						(complex<float> (*)[3]) device_zhelem,
						(complex<float> (*)[3]) device_zgelem,
						(double (*)[4])          device_co,
						cxm[i],
						cym[i],
						czm[i],
						&device_etas[j],
						*zge,
						*zcs,
						*zcp,
						(double (*)[3])device_delta,
						pi,
						*fr,
						device_gi,
						device_ome,
						*npg
						);
				
				/*
                nonsingd(zhelem, 
                          zgelem, 
                          co, 
                          &cxm[i], 
                          &cym[i], 
                          &czm[i],
                          etas[j],
                          zge,
                          zcs,
                          zcp,
                          delta,
                          pi,
                          fr,
                          npg
                         );
			*/
	/*	
				error = cudaMemcpy(zgelemf, device_zgelem, 3*3*sizeof(complex<float>), cudaMemcpyDeviceToHost);
				cuda_assert(error);
				error = cudaMemcpy(zhelemf, device_zhelem, 3*3*sizeof(complex<float>), cudaMemcpyDeviceToHost);
				cuda_assert(error);
*/
				for (int jjj = 0; jjj < 3; ++jjj)
				{   for (int iii = 0; iii < 3; ++iii)
                    {
                        zgp[jjj+jj][iii+ii] = zgelemf[jjj][iii];
                        zhp[jjj+jj][iii+ii] = zhelemf[jjj][iii];
                    }
				}
			}
        }
    }
	printf("ACABOU\n");
}
}
