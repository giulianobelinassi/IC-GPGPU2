#ifndef SHARED_H
#define SHARED_H

#include "cublas_v2.h"

extern "C"{            

__constant__ float delta[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                       
extern float* device_cx;
extern float* device_cy;                                                              
extern float* device_cz;
extern float* device_cxm;
extern float* device_cym;                                                               
extern float* device_czm;
extern float* device_etas;
extern float* device_gi;
extern float* device_ome;
extern int*   device_cone;        
 
void cuda_assert(cudaError_t error);

void cublas_assert(cublasStatus_t error);

int largest_possible_width(size_t sizeof_column_mem, int columns, int* iterations);
}
#endif

