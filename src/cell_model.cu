#include <iostream>
#include <curand_kernel.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#include "cell_model.cuh"
#include "cell_model_params.h"
#include "cuda_util.h"
#include "cell/movement.cuh"

#include <curand_kernel.h>

// TODO: Refactor cuda functions to dedicated files:

// Cuda simulation functions:
__global__ void initialise_grid(GridElement *grid, CellModelParams params);
__global__ void update_cells(GridElement *grid, CellModelParams params);
__global__ void update_environment(GridElement *grid, CellModelParams params);
__global__ void update_interactions(GridElement *grid, CellModelParams params);
__device__ void iterate_cell(GridElement *grid, CellModelParams params, int idx, int x, int y, int z);

// GridElement initialisation functions:
__device__ void initialise_cell(Cell &cell, CellModelParams params);
__device__ void initialise_environement(Environment &environment, CellModelParams params);

// Cell operations:
__device__ void energy_death(GridElement &element, CellModelParams &params);
__device__ void use_energy(GridElement &element, int y, CellModelParams &params);
__device__ void acquire_energy(GridElement &element, int y, CellModelParams &params);

// Util
__device__ int get_idx(int x, int y, int z, CellModelParams &params);

// Environment distributions:
__device__ double light_distribution(int y, int height);
__device__ double co2_distribution(int y, int height);
__device__ double chemical_distribution(int y, int height);

CellModel::CellModel(CellModelParams params) :
		params(params),
		numBlocks(params.cudaParams.numBlocks),
		blockSize(params.cudaParams.blockSize) {
	cudaMallocManaged(&grid, params.gridSize * sizeof(GridElement));
	checkCudaError(cudaPeekAtLastError());
	initialise();
}

void CellModel::writeFrame(int dIdx) {
	int nPixels = params.w * params.h;
	unsigned char* framebuffer = new unsigned char[nPixels * 3];

	for (int i = 0; i < nPixels; i++) {
		if (grid[i + dIdx*nPixels].cell.alive) {
			framebuffer[i*3] = grid[i + dIdx*nPixels].cell.energy;
		} else {
			framebuffer[i*3] = 0;
		}
		framebuffer[i*3 + 1] = 0;
		framebuffer[i*3 + 2] = 0;
	}
	std::cout.write(reinterpret_cast<char *>(framebuffer), nPixels * 3);
}

void CellModel::initialise() {
	initialise_grid<<<numBlocks, blockSize>>>(
		grid,
		params
	);
	checkCudaError(cudaPeekAtLastError());
}

void CellModel::simulate(int nIterations) {
	for (int i = 0; i < nIterations; i++) {
		// Cells:
		update_cells<<<numBlocks, blockSize>>>(grid, params);
		checkCudaError(cudaPeekAtLastError());
		// iterateRandomMovement();

		// Environment
		update_environment<<<numBlocks, blockSize>>>(grid, params);
		checkCudaError(cudaPeekAtLastError());

		// Simulation:
		update_interactions<<<numBlocks, blockSize>>>(grid, params);
	}
}

void CellModel::synchronizeData() {
	checkCudaError(cudaPeekAtLastError());
	cudaDeviceSynchronize();
	checkCudaError(cudaPeekAtLastError());
}

__global__
void initialise_grid(GridElement *grid, CellModelParams params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Stride loop:
	for (int idx = tid; idx < params.gridSize; idx += blockDim.x * gridDim.x) {
		grid[idx] = {};
		curand_init(1000, idx, 0, &grid[idx].randState);

		if (curand_uniform(&grid[idx].randState) < params.initialDensity) {
			Cell &cell = grid[idx].cell;
			initialise_cell(cell, params);
		}

		Environment &environment = grid[idx].environment;
		initialise_environement(environment, params);
	}
}

__device__
void initialise_cell(Cell &cell, CellModelParams params) {
	cell.alive = true;
	cell.energy = 200;
	cell.waste = 50;
}

__device__
void initialise_environement(Environment &environment, CellModelParams params) {

}

__global__
void update_cells(GridElement *grid, CellModelParams params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Stride loop:
	for (int idx = tid; idx < params.gridSize; idx += blockDim.x * gridDim.x) {
		int x = idx % params.w;
		int y = (idx / params.w) % params.h;
		int z = idx / (params.w * params.h);
		iterate_cell(grid, params, idx, x, y, z);
	}
}

__global__
void update_environment(GridElement *grid, CellModelParams params) {
	// int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	// // Stride loop:
	// for (int idx = tid; idx < params.gridSize; idx += blockDim.x * gridDim.x) {
	// 	int x = idx % params.w;
	// 	int y = (idx / params.w) % params.h;
	// 	int z = idx / (params.w * params.h);
	// }
}

__global__
void update_interactions(GridElement *grid, CellModelParams params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Stride loop:
	for (int idx = tid; idx < params.gridSize; idx += blockDim.x * gridDim.x) {
		int x = idx % params.w;
		int y = (idx / params.w) % params.h;
		int z = idx / (params.w * params.h);
		move_cells(grid, grid[idx], x, y, z, params);
	}
}

__device__
void iterate_cell(GridElement *grid, CellModelParams params, int idx, int x, int y, int z) {
	GridElement &element = grid[idx];

	if (element.cell.alive) {
		// Energy:
		energy_death(element, params);
		use_energy(element, y, params);
		acquire_energy(element, y, params);

		// Movement:
		set_move_direction(grid, element, x, y, z, params);
	}
}

__device__
void energy_death(GridElement &element, CellModelParams &params) {
	if (element.cell.energy < params.survivalThreshold) {
		element.cell.alive = false;
		// TODO: Release resources (90%)
	}
}

__device__
void use_energy(GridElement &element, int y, CellModelParams &params) {
	int newEnergy = element.cell.energy - params.energyUsageRate;
	element.cell.energy = newEnergy > 0 ? newEnergy : 0;
}

__device__
void acquire_energy(GridElement &element, int y, CellModelParams &params) {
	int newEnergy = element.cell.energy;

	newEnergy += params.lightEnergyConversionRate * light_distribution(y, params.h);
	newEnergy += params.co2EnergyConversionRate * co2_distribution(y, params.h);

	int maxEnergy = (1 << Cell::nbits_energy) - 1;
	element.cell.energy = newEnergy < maxEnergy ? newEnergy : maxEnergy;
}

__device__
int get_idx(int x, int y, int z, CellModelParams &params) {
	return x + y * params.w + z * params.w * params.h;
}

__device__
double light_distribution(int y, int height) {
	// TODO: height must be greater than 1
	double depth = y/(double)(height - 1);
	return 1 - (depth * depth);
}

__device__
double co2_distribution(int y, int height) {
	double depth = y/(double)(height - 1);
	return 1 - depth;
}

__device__
double chemical_distribution(int y, int height) {
	float depth = y/(float)(height - 1);
	return pow(depth, 3);
}

// --------------------------------
// STATISTICS ---------------------

struct CellOccupied {
	__device__
	int operator()(const GridElement& g) const {
		return g.cell.alive ? 1 : 0;
	}
};

struct CellEnergy {
	__device__
	double operator()(const GridElement& g) const {
		return g.cell.alive ? g.cell.energy : 0.0;
	}
};

int CellModel::numberOfLivingCells() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellOccupied(), 0, thrust::plus<int>());
}

double CellModel::totalEnergy() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellEnergy(), 0.0, thrust::plus<double>());
}


// Light function (quadratic decrease with depth)
// CO2 function (linear decrease with depth)
// Initialise chemical distribution (exponential increase with depth)
// Store energy (probability at energy > threshold, e.g. 150)
// Temperature function (make temperature change energy usage)
// Capture chemicals (costs energy, chem in env. must be 4-way adjacent or coincident)
// Toxins: digestible, indigestible (dies if exceed threshold)
// Digestible toxins: Remove using energy, or export to adjacent cell or environment


// When a GridElement has an interaction with another GridElement, modifications
// must happen in a separate kernel invocation. Modifications to simulation state
// must also happen in separate kernel invocation.