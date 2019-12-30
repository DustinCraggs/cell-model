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

CellModel::CellModel(CellModelParams params) :
		params(params),
		blockSize(params.nThreads, params.nThreads, 1),
		numBlocks((
			params.w + params.nThreads - 1) / params.nThreads,
			(params.h + params.nThreads - 1) / params.nThreads,
			1
		) {
	cudaMallocManaged(&cells, params.w * params.h * sizeof(Cell));
	initialise();
}

void CellModel::printCells() {
	for (int i = 0; i < params.h; i++) {
		for (int j = 0; j < params.w; j++) {
			std::cout << cells[j + i * params.w] << " | ";
		}
		std::cout << std::endl;
	}
}

void CellModel::initialise() {
	initialise_cells<<<numBlocks, blockSize>>>(
		cells,
		params.w,
		params.h,
		params.initial_density
	);
	checkCudaError(cudaPeekAtLastError());
	cudaDeviceSynchronize();
}

void CellModel::simulate(int nIterations) {
	simulate_cells<<<numBlocks, blockSize>>>(cells, params, nIterations);
	checkCudaError(cudaPeekAtLastError());
	cudaDeviceSynchronize();
}

__global__
void initialise_cells(Cell *cells, int w, int h, float density) {
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;

	if (x < w && y < h) {
		int idx = x + y * w;
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
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	if (x < params.w && y < params.h) {
		for (int i = 0; i < nIterations; i++) {
			iterate(cells, params, x, y);
		}
	}
}

__device__
void iterate(Cell *cells, CellModelParams params, int x, int y) {
	int idx = x + y * params.w;
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