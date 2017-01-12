#include <complex>

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
/*  int n1, n2, n3, n4; */
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
/*
		n1 = cone[0][j];
        n2 = cone[1][j];
        n3 = cone[2][j];
        n4 = cone[3][j];
*/
/*
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
*/
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
                    for (int iii = 0; iii < 3; ++iii)
                    {
                        zgp[jjj+jj][iii+ii] = zgelem[jjj][iii] + zgest[jjj+jj][iii+ii];
                        zhp[jjj+jj][iii+ii] = zhelem[jjj][iii] + zhest[jjj+jj][iii+ii];
                    }
            }
            else
            {
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
                    for (int iii = 0; iii < 3; ++iii)
                    {
                        zgp[jjj+jj][iii+ii] = zgelem[jjj][iii];
                        zhp[jjj+jj][iii+ii] = zhelem[jjj][iii];
                    }
            }
        }
    }
}
}
