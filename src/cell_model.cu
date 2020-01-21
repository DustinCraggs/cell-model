#include <iostream>
#include <curand_kernel.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#include "cell_model.cuh"
#include "cell.h"
#include "cuda_util.h"

#include <curand_kernel.h>


__global__ void initialise_cells(Cell *cells, CellModelParams params);
__global__ void simulate_cells(Cell *cells, CellModelParams params, int nIterations);
__device__ void iterate(Cell *cells, CellModelParams params, int idx, int x, int y, int z);
__device__ void energy_death(Cell *cells, int idx, int energy_threshold);
__device__ void use_energy(Cell *cells, int idx, float probability);
__device__ void acquire_energy(Cell *cells, int idx, float probability);

CellModel::CellModel(CellModelParams params) :
		params(params),
		numBlocks(params.cudaParams.numBlocks),
		blockSize(params.cudaParams.blockSize) {
	cudaMallocManaged(&cells, params.nCells * sizeof(Cell));
	initialise();

	std::cout << "Cell size: " << sizeof(Cell) << std::endl;
	std::cout << "Rand state size: " << sizeof(curandState) << std::endl;
	std::cout << "Cells size: " << sizeof(Cell) * params.nCells << std::endl;
}

void CellModel::printCells() {
	for (int i = 0; i < params.d; i++) {
		for (int j = 0; j < params.h; j++) {
			for (int k = 0; k < params.w; k++) {
				std::cout << cells[k + j * params.w + i * params.w * params.h] << " | ";
			}
			std::cout << std::endl;
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
}

void CellModel::simulate(int nIterations) {
	simulate_cells<<<numBlocks, blockSize>>>(cells, params, nIterations);
	checkCudaError(cudaPeekAtLastError());
}

void CellModel::synchronizeData() {
	cudaDeviceSynchronize();
}

struct CellOccupied {
	__device__
	int operator()(const Cell& c) const {
		return c.occupied;
	}
};

struct CellEnergy {
	__device__
	double operator()(const Cell& c) const {
		return c.occupied ? c.energy : 0.0;
	}
};

int CellModel::numberOfLivingCells() {
	thrust::device_ptr<Cell> cellsPtr = thrust::device_pointer_cast(cells);
	return thrust::transform_reduce(cellsPtr, cellsPtr + params.nCells, CellOccupied(), 0, thrust::plus<int>());
}

double CellModel::totalEnergy() {
	thrust::device_ptr<Cell> cellsPtr = thrust::device_pointer_cast(cells);
	return thrust::transform_reduce(cellsPtr, cellsPtr + params.nCells, CellEnergy(), 0.0, thrust::plus<double>());
}

__global__
void initialise_cells(Cell *cells, CellModelParams params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Stride loop:
	for (int idx = tid; idx < params.nCells; idx += blockDim.x * gridDim.x) {
		curand_init(1000, idx, 0, &cells[idx].randState);
		if (curand_uniform(&cells[idx].randState) < params.initialDensity) {
			cells[idx].occupied = true;
			cells[idx].energy = 10;
			cells[idx].waste = 12;
		}
	}
}

__global__
void simulate_cells(Cell *cells, CellModelParams params, int nIterations) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	for (int i = 0; i < nIterations; i++) {
		// Stride loop:
		for (int idx = tid; idx < params.nCells; idx += blockDim.x * gridDim.x) {
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
		use_energy(cells, idx, params.energyLossProb);
		acquire_energy(cells, idx, params.gatherLightEnergyProb);
		energy_death(cells, idx, params.survivalThreshold);
	}
}

__device__
void energy_death(Cell *cells, int idx, int energyThreshold) {
	if (cells[idx].energy < energyThreshold) {
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