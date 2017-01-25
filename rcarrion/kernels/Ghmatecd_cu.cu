#include <complex.h>
#include <cstdio>
#include <cmath>

#define complex _Complex
#define restrict __restrict__

extern "C"{

void gauleg(double* x1, double* x2, double x[], double w[], int* n)
{
	const double eps = 3E-14;
	const double pi  = 3.141592654;
	int m = (*n+1)/2;
	double xm = 0.5*((*x2) + (*x1));
	double xl = 0.5*((*x2) - (*x1));
	int i, j;

	for (i = 0; i < m; ++i)
	{
		double z = cos(pi*((i+1) - 0.25)/((*n) + 0.5));
		double p1, p2, pp, z1, z2;

		do
		{
			p1 = 1.0;
			p2 = 0.0;

			for (j = 0; j < (*n); ++j)
			{
				double p3 = p2;
				p2 = p1;
				p1 = ((2.0*(j) + 1.0)*z*p2 - (j)*p3)/(j+1);
			}
			pp = (*n)*(z*p1 - p2)/(z*z-1.0);
			z1 = z;
			z = z1 - p1/pp;
		} while(fabs(z - z1) > eps);
	
		x[i] = xm - xl*z;
		x[(*n)-i-1] = xm+xl*z;
		w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
		w[(*n)-i-1] = w[i];
	}

}

void solfundif_(double complex zudif[3][3], 
			    double complex ztdif[3][3], 
			    double* cxp, 
			    double* cyp, 
			    double* czp, 
			    double* cxg, 
			    double* cyg, 
			    double* czg, 
			    double rn[3], 
			    double complex* zge,
			    double complex* zcs, 
			    double complex* zcp,
				double* c1,
				double* c2,
				double* c3,
				double* c4,
			    double delta[3][3], 
			    double* pi, 
			    double* fr
			 );

void solfundif (double complex zudif[3][3], 
			    double complex ztdif[3][3], 
			    double* cxp, 
			    double* cyp, 
			    double* czp, 
			    double* cxg, 
			    double* cyg, 
			    double* czg, 
			    double rn[3], 
			    double complex* zge,
			    double complex* zcs, 
			    double complex* zcp,
				double* c1,
				double* c2,
				double* c3,
				double* c4,
			    double delta[3][3], 
			    double* pi, 
			    double* fr
			 )
{
	
    double r1, r2, r3, r, drn, rd[3];
    double complex zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                   zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;
    
    int i, j;

    r1    = *cxg - *cxp;
    r2    = *cyg - *cyp;
    r3    = *czg - *czp;
    r     = sqrt(r1*r1 + r2*r2 + r3*r3);
    drn   = (r1*rn[0] + r2*rn[1] + r3*rn[2])/r;
    rd[0] = r1/r;
    rd[1] = r2/r;
    rd[2] = r3/r;

    zwi = 0 + (*fr) * (_Complex_I);
    
    zc0 = 1.0/(4*(*pi)*(*zge));
    zc1 = ((*zcp)/(*zcs))*((*zcp)/(*zcs));
    zc2 = ((*zcs)/(*zcp))*((*zcs)/(*zcp));
    zkp = -zwi/(*zcp);
    zks = -zwi/(*zcs);
    zzp = zkp*r;
    zzs = zks*r;
    zezp= cexp(zzp);
    zezs= cexp(zzs);
    zp2 = zzp*zzp;
    zs2 = zzs*zzs;

    zfhi    = (1 + 1./zs2 - 1./zzs)*zezs/r - zc2*(1./zp2 - 1./zzp)*zezp/r;
    zcappa  = (1 + 3./zs2 - 3./zzs)*zezs/r - zc2*(1 + 3./zp2 - 3./zzp)*zezp/r;
    zfhidr  = (-2+ zzs + 3./zzs - 3./zs2)*zezs/(r*r) - zc2*(-1 + 3./zzp - 3./zp2)*zezp/(r*r);
    zcappadr= (zzs - 4. + 9./zzs - 9./zs2)*zezs/(r*r) - zc2*(zzp - 4. + 9./zzp - 9./zp2)*zezp/(r*r);

    zaa = zfhidr-zcappa/r;
    zbb = 4.*zcappa/r - 2.*zcappadr;
    zcc = (zc1-2.)*(zaa + 0.5*zbb-3.0*zcappa/r)-2.0*zcappa/r;

    for (j = 0; j < 3; ++j)
    {   for (i = 0; i < 3; ++i)
        {
            zudif[j][i] = zc0*(zfhi*delta[j][i] - zcappa*rd[j]*rd[i]) - ((*c1)/r)*((*c2)*delta[j][i] + rd[i]*rd[j]);
//        	ztdif[j][i] = (1.0/(4.0*(*pi)))*((zaa*(drn*delta[j][i]+rd[j]*rn[i])) + rd[i]*rd[j]*drn*zbb + rd[i]*rn[j]*zcc) + ((*c3)/(r*r))*(rdn*((*c4)*delta[j][i] + 3.0*rd[i]*rd[j] + (*c4)*(rd[j]*rn[i] - rd[i]*rn[j]));	
        	ztdif[j][i] = (1.0/(4.0*(*pi)))*((zaa*(drn*delta[j][i]+rd[j]*rn[i])) +
					rd[i]*rd[j]*drn*zbb + rd[i]*rn[j]*zcc) + 
					((*c3)/(r*r))*(drn*((*c4)*delta[j][i] + 
					3.0*rd[i]*rd[j]) + 
					(*c4)*(rd[j]*rn[i] - rd[i]*rn[j]));
		}
    }
}

void sing_de_(double complex zhelem[3][3],
              double complex zgelem[3][3],
              double co[3][4],
              double* cxm,
              double* cym,
              double* czm,
              double eta[3],
              double complex* zge,
              double complex* zcs,
              double complex* zcp,
              double* c1,
              double* c2,
              double* c3,
              double* c4,
              double delta[3][3],
              double* pi,
              double* fr,
              int* npg
              );



void sing_de( double complex zhdifel[3][3],
              double complex zgdifel[3][3],
              double co[3][4],
              double* cxp,
              double* cyp,
              double* czp,
              double eta[3],
              double complex* zge,
              double complex* zcs,
              double complex* zcp,
              double* c1,
              double* c2,
              double* c3,
              double* c4,
              double delta[3][3],
              double* pi,
              double* fr,
              int* npg
              )
{
    double complex zudif[3][3];
    double complex ztdif[3][3];
    double gi[*npg], ome[*npg], p[4][2], xj[3][2], f[4];
	double g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
    double cxg, cyg, czg;
    double j1, j2, j3;
	double minus1 = -1;
	double plus1  = +1;
    int i, j, k, ig, jg;

    for (j = 0; j < 3; ++j)
        for (i = 0; i < 3; ++i)
        {
            zhdifel[j][i] = 0;
            zgdifel[j][i] = 0;
        }

    gauleg(&minus1, &plus1, gi , ome , npg);

    for (jg = 0; jg < *npg; ++jg)
    {
        g2 = gi[jg];
        p2 = ome[jg];
        sp = 1 + g2;
        sm = 1 - g2;
        p[0][0] = -0.25*sm;
        p[1][0] =  0.25*sm;
        p[2][0] =  0.25*sp;
        p[3][0] = -0.25*sp;

        for (ig = 0; ig < *npg; ++ig)
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

            for (i = 0; i < 2; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    temp = 0;
                    for (k = 0; k < 4; ++k)
                        temp = temp + p[k][i]*co[j][k];
                    xj[j][i] = temp;
                }
            }

            j1 = xj[1][0]*xj[2][1]-xj[1][1]*xj[2][0];
            j2 = xj[0][1]*xj[2][0]-xj[0][0]*xj[2][1];
            j3 = xj[0][0]*xj[1][1]-xj[0][1]*xj[1][0];

            det = sqrt(j1*j1 + j2*j2 + j3*j3);

            if (det < 1e-5)
            {
                printf("JACOBIANO NULO OU NEGATIVO!!!\n");
            }

			cxg = 0;
			cyg = 0;
			czg = 0;
			for (i = 0; i < 4; ++i)
			{
				cxg = cxg + co[0][i]*f[i];
				cyg = cyg + co[1][i]*f[i];
				czg = czg + co[2][i]*f[i];
			}

			solfundif_(zudif, ztdif, cxp, cyp, czp, &cxg, &cyg, &czg, eta, zge, zcs, zcp, c1, c2, c3, c4, delta, pi, fr);

			p12 = p1*p2*det;
			for (j = 0; j < 3; ++j)
			{
				for (i = 0; i < 3; ++i)
				{
					zhdifel[j][i] = zhdifel[j][i] + ztdif[j][i]*p12;
					zgdifel[j][i] = zgdifel[j][i] + zudif[j][i]*p12;
				}
			}
        }
    }
}

void nonsingd_(double complex zhelem[3][3],
               double complex zgelem[3][3],
              double co[3][4],
              double* cxm,
              double* cym,
              double* czm,
              double eta[3],
              double complex* zge,
              double complex* zcs,
              double complex* zcp,
              double delta[3][3],
              double* pi,
              double* fr,
              int* npg
              );



void gauleg_(double* x1, double* x2, double x[], double w[], int* n);


void solfund_(double complex zu[3][3], 
			  double complex zt[3][3], 
			  double* cxp, 
			  double* cyp, 
			  double* czp, 
			  double* cxg, 
			  double* cyg, 
			  double* czg, 
			  double rn[3], 
			  double complex* zge,
			  double complex* zcs, 
			  double complex* zcp, 
			  double delta[3][3], 
			  double* pi, 
			  double* fr
			 );

void solfund(double complex zu[3][3], 
			  double complex zt[3][3], 
			  double* cxp, 
			  double* cyp, 
			  double* czp, 
			  double* cxg, 
			  double* cyg, 
			  double* czg, 
			  double rn[3], 
			  double complex* zge,
			  double complex* zcs, 
			  double complex* zcp, 
			  double delta[3][3], 
			  double* pi, 
			  double* fr
			 )
{
    double r1, r2, r3, r, drn, rd[3];
    double complex zwi, zc0, zc1, zc2, zkp, zks, zzp, zzs, zezp, zezs, 
                   zp2, zs2, zfhi, zcappa, zfhidr, zcappadr, zaa, zbb, zcc;
    
    int i, j;

    r1    = *cxg - *cxp;
    r2    = *cyg - *cyp;
    r3    = *czg - *czp;
    r     = sqrt(r1*r1 + r2*r2 + r3*r3);
    drn   = (r1*rn[0] + r2*rn[1] + r3*rn[2])/r;
    rd[0] = r1/r;
    rd[1] = r2/r;
    rd[2] = r3/r;

    zwi = 0 + (*fr) * (_Complex_I);
    
    zc0 = 1.0/(4*(*pi)*(*zge));
    zc1 = ((*zcp)/(*zcs))*((*zcp)/(*zcs));
    zc2 = ((*zcs)/(*zcp))*((*zcs)/(*zcp));
    zkp = -zwi/(*zcp);
    zks = -zwi/(*zcs);
    zzp = zkp*r;
    zzs = zks*r;
    zezp= cexp(zzp);
    zezs= cexp(zzs);
    zp2 = zzp*zzp;
    zs2 = zzs*zzs;

    zfhi    = (1 + 1./zs2 - 1./zzs)*zezs/r - zc2*(1./zp2 - 1./zzp)*zezp/r;
    zcappa  = (1 + 3./zs2 - 3./zzs)*zezs/r - zc2*(1 + 3./zp2 - 3./zzp)*zezp/r;
    zfhidr  = (-2+ zzs + 3./zzs - 3./zs2)*zezs/(r*r) - zc2*(-1 + 3./zzp - 3./zp2)*zezp/(r*r);
    zcappadr= (zzs - 4. + 9./zzs - 9./zs2)*zezs/(r*r) - zc2*(zzp - 4. + 9./zzp - 9./zp2)*zezp/(r*r);

    zaa = zfhidr-zcappa/r;
    zbb = 4.*zcappa/r - 2.*zcappadr;
    zcc = (zc1-2.)*(zaa + 0.5*zbb-3.0*zcappa/r)-2.0*zcappa/r;

    for (j = 0; j < 3; ++j)
    {   for (i = 0; i < 3; ++i)
        {
            zu[j][i] = zc0*(zfhi*delta[j][i] - zcappa*rd[j]*rd[i]);
        	zt[j][i] = (1.0/(4.0*(*pi)))*((zaa*(drn*delta[j][i]+rd[j]*rn[i])) + rd[i]*rd[j]*drn*zbb + rd[i]*rn[j]*zcc);
		}
    }
}



void nonsingd(double complex zhelem[3][3],
              double complex zgelem[3][3],
              double co[3][4],
              double* cxp,
              double* cyp,
              double* czp,
              double eta[3],
              double complex* zge,
              double complex* zcs,
              double complex* zcp,
              double delta[3][3],
              double* pi,
              double* fr,
              int* npg
              )
{
    double complex zu[3][3];
    double complex zt[3][3];
    double gi[*npg], ome[*npg], p[4][2], xj[3][2], f[4];
	double g1, g2, p1, p2, p12, sp, sm, rp, rm, temp, det;
    double cxg, cyg, czg;
    double j1, j2, j3;
	double minus1 = -1;
	double plus1  = +1;
    int i, j, k, ig, jg;

    for (j = 0; j < 3; ++j)
        for (i = 0; i < 3; ++i)
        {
            zhelem[j][i] = 0;
            zgelem[j][i] = 0;
        }

    gauleg(&minus1, &plus1, gi , ome , npg);

    for (jg = 0; jg < *npg; ++jg)
    {
        g2 = gi[jg];
        p2 = ome[jg];
        sp = 1 + g2;
        sm = 1 - g2;
        p[0][0] = -0.25*sm;
        p[1][0] =  0.25*sm;
        p[2][0] =  0.25*sp;
        p[3][0] = -0.25*sp;

        for (ig = 0; ig < *npg; ++ig)
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

            for (i = 0; i < 2; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    temp = 0;
                    for (k = 0; k < 4; ++k)
                        temp = temp + p[k][i]*co[j][k];
                    xj[j][i] = temp;
                }
            }

            j1 = xj[1][0]*xj[2][1]-xj[1][1]*xj[2][0];
            j2 = xj[0][1]*xj[2][0]-xj[0][0]*xj[2][1];
            j3 = xj[0][0]*xj[1][1]-xj[0][1]*xj[1][0];

            det = sqrt(j1*j1 + j2*j2 + j3*j3);

            if (det < 1e-5)
            {
                printf("JACOBIANO NULO OU NEGATIVO!!!\n");
            }

			cxg = 0;
			cyg = 0;
			czg = 0;
			for (i = 0; i < 4; ++i)
			{
				cxg = cxg + co[0][i]*f[i];
				cyg = cyg + co[1][i]*f[i];
				czg = czg + co[2][i]*f[i];
			}

			solfund(zu, zt, cxp, cyp, czp, &cxg, &cyg, &czg, eta, zge, zcs, zcp, delta, pi, fr);

			p12 = p1*p2*det;
			for (j = 0; j < 3; ++j)
			{
				for (i = 0; i < 3; ++i)
				{
					zhelem[j][i] = zhelem[j][i] + zt[j][i]*p12;
					zgelem[j][i] = zgelem[j][i] + zu[j][i]*p12;
				}
			}
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
                    double complex* zge,
                    double complex* zcs,
                    double complex* zcp,
                    double* c1,
                    double* c2,
                    double* c3,
                    double* c4,
                    double delta[3][3],
                    double* pi,
                    double* fr,
                    double* zhest_,
					double* zgest_,
					double complex* zgp_,
                    double complex* zhp_,
					int* status
                   )
{
    double co[3][4];
    int n1, n2, n3, n4;
    int i, j, ii, jj;
    double complex zgelem[3][3];
	double complex zhelem[3][3];

	/*Cast os parâmetros de volta para o tipo original*/
	int (*cone)[*ne]           = (int (*)[*ne]) cone_;
	double (*zhest)[*nx]       = (double (*)[*nx]) zhest_;
	double (*zgest)[*nx]       = (double (*)[*nx]) zgest_;
	double complex (*zgp)[*nx] = (double complex (*)[*nx]) zgp_;
	double complex (*zhp)[*nx] = (double complex (*)[*nx]) zhp_;
/*	
	double co_cpy[3][4];
	int cone_cpy[4][*ne];
	double cxm_cpy[*ne], cym_cpy[*ne], czm_cpy[*ne];
	double eta_cpy[3];
	double complex zge_cpy, zcs_cpy, zcp_cpy;
	double delta_cpy[3][3];
	double pi_cpy, fr_cpy;
	int npg_cpy;
*/

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

        jj = 3*j;
        for (i = 0; i < *nbe; ++i)
        {
            ii = 3*i;

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
                         pi,
                         fr,
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
/*				
				for (int y = 0; y < 3; ++y)
				{
					for (int x = 0; x < 4; ++x)
					{
						co_cpy[y][x] = co[y][x];
					}
				}

				for (int x = 0; x < *ne; ++x)
				{
					cxm_cpy[x] = cxm[x];
					cym_cpy[x] = cym[x];
					czm_cpy[x] = czm[x];
				}

				for (int x = 0; x < 3; ++x)
				{
					eta_cpy[x] = etas[j][x];
				}
	
				zge_cpy = *zge;
				zcs_cpy = *zcs;
				zcp_cpy = *zcp;
				
				for (int y = 0; y < 3; ++y)
					for (int x = 0; x < 3; ++x)
						delta_cpy[y][x] = delta[y][x];

				pi_cpy = *pi;
				fr_cpy = *fr;
				npg_cpy = *npg;
*/
                nonsingd_(zhelem, 
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
			
				for (int jjj = 0; jjj < 3; ++jjj)
				{   for (int iii = 0; iii < 3; ++iii)
                    {
                        zgp[jjj+jj][iii+ii] = zgelem[jjj][iii];
                        zhp[jjj+jj][iii+ii] = zhelem[jjj][iii];
                    }
				}

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
/*
				for (int y = 0; y < 3; ++y)
				{
					for (int x = 0; x < 4; ++x)
					{
						if (co[y][x] != co_cpy[y][x])
							printf("CO ALTERADO!, %d, %d\n", i, j);
					}
				}

				for (int x = 0; x < *ne; ++x)
				{
					if (cxm_cpy[x] != cxm[x])
						printf("CXM ALTERADO!, %d, %d\n", i, j);
					if (cym_cpy[x] != cym[x])
						printf("CYM ALTERADO!, %d, %d\n", i, j);
					if (czm_cpy[x] != czm[x])
						printf("CZM ALTERADO!, %d, %d\n", i, j);
				}

				for (int x = 0; x < 3; ++x)
				{
					if (eta_cpy[x] != etas[j][x])
						printf("ETA ALTERADO!");
				}


				if (zge_cpy != *zge)
					printf("ZGE ALTERADO!\n");
				if (zcs_cpy != *zcs)
					printf("ZCS ALTERADO!\n");
				if (zcp_cpy != *zcp)
					printf("ZCP ALTERADO!\n");

				if (fr_cpy != *fr)
					printf("FR ALTERADO!\n");
				if (pi_cpy != *pi)
					printf("PI ALTERADO!\n");
				if (npg_cpy != *npg)
					printf("NPG ALTERADO!\n");
*/				
                for (int jjj = 0; jjj < 3; ++jjj)
				{   for (int iii = 0; iii < 3; ++iii)
                    {
                        if (zgp[jjj+jj][iii+ii] != zgelem[jjj][iii])
							printf("ZG i: %d, j: %d\n", i, j);
                        if (zhp[jjj+jj][iii+ii] != zhelem[jjj][iii])
							printf("ZH i: %d, j: %d\n", i, j);
                    }
				}

			}
        }
    }
}
}
