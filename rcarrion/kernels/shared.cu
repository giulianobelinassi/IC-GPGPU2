#include "shared.h"
#include <cstdio>

extern "C"{

float* device_cx;
float* device_cy;
float* device_cz;
float* device_cxm;
float* device_cym;
float* device_czm;
float* device_etas;
float* device_gi;
float* device_ome;
int*   device_cone;

void cuda_assert(cudaError_t error)
{
    if (error != cudaSuccess)
    {   
        fputs(cudaGetErrorString(cudaGetLastError()), stderr);
        putc('\n', stderr);
        exit(1);
    }   
}

void send_shared_data_to_gpu_(
        float cx[],
        float cy[],
        float cz[],
        float cxm[],
        float cym[],
        float czm[],
        float etas[],
        float gi[],
        float ome[],
        int cone[],
        int* np, 
        int* npg,
        int* n,
        int* nbe 
        )
{
    cudaError_t error;

    /*Aloque memória para os vetores na GPU*/
    error = cudaMalloc(&device_cone, 4*(*n)*sizeof(int));
    cuda_assert(error);

    error = cudaMalloc(&device_cx, (*np)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_cy, (*np)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_cz, (*np)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_cxm, (*n)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_cym, (*n)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_czm, (*n)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_gi, (*npg)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_ome, (*npg)*sizeof(float));
    cuda_assert(error);

    error = cudaMalloc(&device_etas, (*n)*3*sizeof(float));
    cuda_assert(error);

    /*mova os dados para lá*/
    error = cudaMemcpy(device_cone, cone, 4*(*n)*sizeof(int), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cx, cx, (*np)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cy, cy, (*np)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cz, cz, (*np)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cxm, cxm, (*n)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_cym, cym, (*n)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_czm, czm, (*n)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_gi, gi, (*npg)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_ome, ome, (*npg)*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

    error = cudaMemcpy(device_etas, etas, (*n)*3*sizeof(float), cudaMemcpyHostToDevice);
    cuda_assert(error);

}

void deallocate_shared_gpu_data_()
{
    cudaError_t error;

    error = cudaFree(device_cone);
    cuda_assert(error);
    error = cudaFree(device_gi);
    cuda_assert(error);
    error = cudaFree(device_ome);
    cuda_assert(error);
    error = cudaFree(device_etas);
    cuda_assert(error);
    error = cudaFree(device_cx);
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