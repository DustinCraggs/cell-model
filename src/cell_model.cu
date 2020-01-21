#include <iostream>
#include <curand_kernel.h>

#include "cell_model.cuh"
#include "cell.h"
#include "cuda_util.h"

__global__ void initialise_cells(Cell *cells, int w, int h, float density);
__global__ void simulate_cells(Cell *cells, CellModelParams params, int nIterations);
__device__ void iterate(Cell *cells, CellModelParams params, int x, int y);
__device__ void energy_death(Cell *cells, int idx, int energy_threshold);
__device__ void use_energy(Cell *cells, int idx, float probability);
__device__ void acquire_energy(Cell *cells, int idx, float probability);

CellModel::CellModel(CellModelParams params) : params(params) {
	int nCells = params.w * params.h * params.d;
	cudaMallocManaged(&cells, nCells * sizeof(Cell));
	initialise();
}

void CellModel::printCells() {
	for (int i = 0; i < params.h; i++) {
		for (int j = 0; j < params.w; j++) {
			for (int k = 0;  < params.d; k++) {
				std::cout << cells[j + i * params.w] << " | ";
			}
		}
		std::cout << std::endl;
	}
}

void CellModel::initialise() {
	initialise_cells<<<numBlocks, blockSize>>>(
		cells,
		params
	);
	checkCudaError(cudaPeekAtLastError());
	cudaDeviceSynchronize();
}

void CellModel::simulate(int nIterations) {
	simulate_cells<<<numBlocks, blockSize>>>(cells, params, nIterations);
	checkCudaError(cudaPeekAtLastError());
	cudaDeviceSynchronize();
}

double* CellModel::getStatistics() {

}

__global__
void initialise_cells(Cell *cells, CellModelParams params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int n = params.w * params.h * params.d;

	// Stride loop:
	for (int idx = tid; tid < n; tid += blockDim.x * gridDim.x) {
		curand_init(1000, idx, 0, &cells[idx].randState);
		if (curand_uniform(&cells[idx].randState) < density) {
			cells[idx].occupied = true;
			cells[idx].energy = 10;
			cells[idx].waste = 12;
		}
	}
}

__global__
void simulate_cells(Cell *cells, CellModelParams params, int nIterations) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int n = params.w * params.h * params.d;

	for (int i = 0; i < nIterations; i++) {
		// Stride loop:
		for (int idx = tid; tid < n; tid += blockDim.x * gridDim.x) {
			int x = idx % params.w;
			int y = (idx / params.w) % params.h;
			int z = idx / (params.w * params.h);
			iterate(cells, params, idx, x, y, z);
		}
	}
}

__device__
void iterate(Cell *cells, CellModelParams params, int idx, int x, int y, int z) {
	if (cells[idx].occupied) {
		use_energy(cells, idx, params.energy_loss_p);
		acquire_energy(cells, idx, params.gather_light_energy_p);
		energy_death(cells, idx, params.survival_threshold);
	}
}

__device__
void energy_death(Cell *cells, int idx, int energy_threshold) {
	if (cells[idx].energy < energy_threshold) {
		cells[idx].occupied = false;
	}
}

__device__
void use_energy(Cell *cells, int idx, float probability) {
	if (curand_uniform(&cells[idx].randState) < probability) {
		cells[idx].energy -= 1;
	}
}

__device__
void acquire_energy(Cell *cells, int idx, float probability) {
	if (cells[idx].energy < pow((int)cells[idx].nbits_energy, 2)) {
		if (curand_uniform(&cells[idx].randState) < probability) {
			cells[idx].energy += 1;
		}
	}
}