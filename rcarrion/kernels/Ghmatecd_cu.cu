#include <cstdio>
#include <cmath>
#include "cuda_complex.hpp"

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
	int i = blockIdx.y;
	int j = blockIdx.x;
	int index = 3*blockIdx.y + 3*blockIdx.x*nx + threadIdx.x*nx + threadIdx.y;
	
	double p[4][2], xj[3][2], f[4];
	double g0, g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
    double cxg, cyg, czg, cxp, cyp, czp;
    double j1, j2, j3;
	double co[3][4];
    double r1, r2, r3, r, drn, rd[3];
    complex<double> zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, zr, zr2, 
                    zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;

	
	const double pi  = 3.141592654;
	
	double temp1, temp2, temp3, temp4;

	int n1 = cone[ne*0 + j];
	int n2 = cone[ne*1 + j];
	int n3 = cone[ne*2 + j];
	int n4 = cone[ne*3 + j];

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


            if (det < 1e-5)
            {
                //algo deveria acontecer
            }


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
	double* device_zhest;
	double* device_zgest;
	complex<double>* device_zh;
	complex<double>* device_zg;
	double* device_cx;
	double* device_cy;
	double* device_cz;
	double* device_cxm;
	double* device_cym;
	double* device_czm;
	double* device_etas;
	double* device_delta;
	int* device_cone;
	double* device_gi;
	double* device_ome;

	int i, j;


	complex<double> (*zgp)[*nx] = (complex<double> (*)[*nx]) zgp_;
	complex<double> (*zhp)[*nx] = (complex<double> (*)[*nx]) zhp_;
	
	
	error = cudaMalloc(&device_gi, (*npg)*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_ome, (*npg)*sizeof(double));
	cuda_assert(error);
	
	error = cudaMalloc(&device_zh, (*nx)*(*nx)*sizeof(complex<double>));
	cuda_assert(error);

	error = cudaMalloc(&device_zg, (*nx)*(*nx)*sizeof(complex<double>));
	cuda_assert(error);
	
	error = cudaMalloc(&device_zhest, (*nx)*(*nx)*sizeof(complex<double>));
	cuda_assert(error);

	error = cudaMalloc(&device_zgest, (*nx)*(*nx)*sizeof(complex<double>));
	cuda_assert(error);

	error = cudaMalloc(&device_etas, (*nx)*3*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_cx, (*ncox)*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_cy, (*ncox)*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_cz, (*ncox)*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_cxm, (*ne)*sizeof(double));	
	cuda_assert(error);

	error = cudaMalloc(&device_cym, (*ne)*sizeof(double));	
	cuda_assert(error);
	
	error = cudaMalloc(&device_czm, (*ne)*sizeof(double));	
	cuda_assert(error);

	error = cudaMalloc(&device_delta, 3*3*sizeof(double));
	cuda_assert(error);

	error = cudaMalloc(&device_cone, 4*(*ne)*sizeof(int));
	cuda_assert(error);

	error = cudaMemcpy(device_gi, gi, (*npg)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_ome, ome, (*npg)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_etas, etas, (*n)*3, cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_cx, cx, (*ncox)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_cy, cy, (*ncox)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_cz, cz, (*ncox)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_cxm, cxm, (*ne)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_cym, cym, (*ne)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_czm, czm, (*ne)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);
	
	error = cudaMemcpy(device_delta, delta, 3*3*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_cone, cone_, (*ne)*4*sizeof(int), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_zhest, zhest_, (*nx)*(*nx)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);

	error = cudaMemcpy(device_zgest, zgest_, (*nx)*(*nx)*sizeof(double), cudaMemcpyHostToDevice);
	cuda_assert(error);


	ghmatecd_kernel<<<numBlocks, threadsPerBlock>>>(*ne,
			                                        *nbe,
													*nx,
													*npg,
													*ncox,
													*n,
													device_zh,
													device_zg,
													(complex<double>) *zge,
													(complex<double>) *zcs,
													(complex<double>) *zcp,
													device_zhest,
													device_zgest,
													device_cx,
													device_cy,
													device_cz,
													device_cxm,
													device_cym,
													device_czm,
													device_gi,
													device_ome,
													(double (*)[3]) device_delta,
													(double (*)[3]) device_etas,
													*fr,
													device_cone
												   );


//	error = cudaMemcpy(zhp_, device_zh, (*nx)*(*nx)*sizeof(complex<double>), cudaMemcpyDeviceToHost);
//	cuda_assert(error);

	error = cudaMemcpy(zgp_, device_zg, (*nx)*(*nx)*sizeof(complex<double>), cudaMemcpyDeviceToHost);
	cuda_assert(error);

	for (j = 0; j < *n; ++j)
	{
		for (i = 0; i < *nbe; ++i)
		{
			printf("%d %d  =   %f\n", i, j, zgp[j][i]);
		}
	}

}
}
