#pragma once
#include <stdio.h>

// https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define checkCudaError(result) { gpuAssert((result), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
   if (code != cudaSuccess) {
      fprintf(stderr, "Cuda error: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}