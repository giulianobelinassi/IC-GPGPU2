#include <complex>
#include <cstdio>

#define complex _Complex
#define restrict __restrict__

extern "C"{

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

void gauleg_(double* , double*, double*, double*, int*);

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

    gauleg_(&minus1, &plus1, gi, ome, npg);

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

			double ccxp = *cxp;
			double ccyp = *cyp;
			double cczp = *czp;

			solfund_(zu, zt, &ccxp, &ccyp, &cczp, &cxg, &cyg, &czg, eta, zge, zcs, zcp, delta, pi, fr);

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
                for (int jjj = 0; jjj < 3; ++jjj)
				{   for (int iii = 0; iii < 3; ++iii)
                    {
                        zgp[jjj+jj][iii+ii] = zgelem[jjj][iii];
                        zhp[jjj+jj][iii+ii] = zhelem[jjj][iii];
                    }
				}
			}
        }
    }
}
}
