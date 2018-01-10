#ifndef SHARED_H
#define SHARED_H

#include "cublas_v2.h"

extern "C"{            

__constant__ FREAL delta[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                       
extern FREAL* device_cx;
extern FREAL* device_cy;                                                              
extern FREAL* device_cz;
extern FREAL* device_cxm;
extern FREAL* device_cym;                                                               
extern FREAL* device_czm;
extern FREAL* device_etas;
extern FREAL* device_gi;
extern FREAL* device_ome;
extern FREAL* device_hestdiag;
extern FREAL* device_gestdiag;
extern int*   device_cone;        
 
void cuda_assert(cudaError_t error);

void cublas_assert(cublasStatus_t error);

int largest_possible_width(size_t sizeof_column_mem, int columns, int* iterations);
}
#endif

