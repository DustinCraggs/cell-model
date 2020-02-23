#include <iostream>
#include <algorithm>
#include <curand_kernel.h>

#include "cell_model.cuh"
#include "cell_model_params.h"
#include "util/cuda_util.h"
#include "operation/movement.cuh"

// TODO: Refactor cuda functions to dedicated files:

// Cuda simulation functions:
__global__ void initialise_grid(GridElement *grid, CellModelParams params);
__global__ void update_cells(GridElement *grid, CellModelParams params);
__global__ void update_environment(GridElement *grid, CellModelParams params);
__global__ void update_interactions(GridElement *grid, CellModelParams params);
__device__ void iterate_cell(GridElement *grid, CellModelParams params, int idx, int x, int y, int z);

// GridElement initialisation functions:
__device__ void initialise_cell(Cell &cell, CellModelParams params);
__device__ void initialise_environement(GridElement &element, int idx, CellModelParams params);

// Cell operations:
__device__ void check_death(GridElement &element, CellModelParams &params);
__device__ void use_energy(GridElement &element, int y, CellModelParams &params);
__device__ void use_chem(GridElement &element, int y, CellModelParams &params);
__device__ void acquire_energy(GridElement &element, int y, CellModelParams &params);
__device__ void acquire_chem(GridElement &element, int x, int y, int z, CellModelParams &params);
__device__ void create_d_toxin(GridElement &element, CellModelParams &params);
__device__ void digest_d_toxin(GridElement &element, CellModelParams &params);

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

GridElement* CellModel::getHostGrid() {
	synchronizeData();
	return grid;
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

		// Environment
		update_environment<<<numBlocks, blockSize>>>(grid, params);
		checkCudaError(cudaPeekAtLastError());

		// Simulation:
		update_interactions<<<numBlocks, blockSize>>>(grid, params);
		checkCudaError(cudaPeekAtLastError());
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

		if (curand_uniform(&grid[idx].randState) < params.initialCellDensity) {
			Cell &cell = grid[idx].cell;
			initialise_cell(cell, params);
		}

		GridElement &element = grid[idx];
		initialise_environement(element, idx, params);
	}
}

__device__
void initialise_cell(Cell &cell, CellModelParams params) {
	cell.alive = true;
	cell.energy = 230;
	cell.chem = 230;
	cell.dToxin = 10;
	cell.ndToxin = 10;
}

__device__
void initialise_environement(GridElement &element, int idx, CellModelParams params) {
	if (curand_uniform(&element.randState) < params.initialChemDensity) {
		int y = (idx / params.w) % params.h;
		element.environment.chem = chemical_distribution(y, params.h) * params.initialChemMax;
	}
	if (curand_uniform(&element.randState) < params.initialNdToxinDensity) {
		element.environment.ndToxin = curand_uniform(&element.randState) * params.initialNdToxinMax;
	}
	element.environment.dToxin = 0;
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

}

__global__
void update_interactions(GridElement *grid, CellModelParams params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Stride loop:
	for (int idx = tid; idx < params.gridSize; idx += blockDim.x * gridDim.x) {
		int x = idx % params.w;
		int y = (idx / params.w) % params.h;
		int z = idx / (params.w * params.h);
		if (grid[idx].cell.alive) {
			move_cells(grid, grid[idx], x, y, z, params);
		}
	}
}

__device__
void iterate_cell(GridElement *grid, CellModelParams params, int idx, int x, int y, int z) {
	GridElement &element = grid[idx];

	if (element.cell.alive) {
		check_death(element, params);
		use_energy(element, y, params);
		use_chem(element, y, params);
		acquire_energy(element, y, params);
		acquire_chem(element, x, y, z, params);
		create_d_toxin(element, params);
		digest_d_toxin(element, params);

		// Movement:
		set_move_direction(grid, element, x, y, z, params);
	}
}

__device__
void check_death(GridElement &element, CellModelParams &params) {
	if (element.cell.energy < params.energySurvivalThreshold
			|| element.cell.chem < params.chemSurvivalThreshold
			|| element.cell.dToxin >= params.dToxinDeathThreshold) {
		element.cell.alive = false;
		// Release 90% of resources to env.:
		int maxChem = (1 << Environment::nbits_chem) - 1;
		int newChem = 0.9 * element.cell.chem + element.environment.chem;
		element.environment.chem = newChem <= maxChem ? newChem : maxChem;

		int maxDToxin = (1 << Environment::nbits_d_toxin) - 1;
		int newDToxin = 0.9 * element.cell.dToxin + element.environment.dToxin;
		element.environment.dToxin = newDToxin <= maxDToxin ? newDToxin : maxDToxin;

		int maxNDToxin = (1 << Environment::nbits_nd_toxin) - 1;
		int newNDToxin = 0.9 * element.cell.ndToxin + element.environment.ndToxin;
		element.environment.ndToxin = newNDToxin <= maxNDToxin ? newNDToxin : maxNDToxin;

		// TODO: Waste, energy
	}
}

__device__
void use_energy(GridElement &element, int y, CellModelParams &params) {
	int newEnergy = element.cell.energy - params.energyUsageRate;
	element.cell.energy = newEnergy > 0 ? newEnergy : 0;
}

__device__
void use_chem(GridElement &element, int y, CellModelParams &params) {
	int newChem = element.cell.chem - params.chemUsageRate;
	element.cell.chem = newChem > 0 ? newChem : 0;
}

__device__
void acquire_energy(GridElement &element, int y, CellModelParams &params) {
	int newEnergy = element.cell.energy;

	newEnergy += params.lightEnergyConversionRate * light_distribution(y, params.h);
	newEnergy += params.co2EnergyConversionRate * co2_distribution(y, params.h);

	int maxEnergy = (1 << Cell::nbits_energy) - 1;
	element.cell.energy = newEnergy <= maxEnergy ? newEnergy : maxEnergy;
}

__device__
void acquire_chem(GridElement &element, int x, int y, int z, CellModelParams &params) {
	int quantity = min(params.chemAcquisitionRate, element.environment.chem);
	
	int maxChem = (1 << Cell::nbits_chem) - 1;
	quantity = min(quantity, maxChem - element.cell.chem);

	element.cell.chem += quantity;
	element.environment.chem -= quantity;
}

__device__
void digest_d_toxin(GridElement &element, CellModelParams &params) {
	if (element.cell.dToxin >= params.dToxinDigestionThreshold) {
		int newDToxin = element.cell.dToxin - params.digestibleToxinDigestionRate;
		element.cell.dToxin = newDToxin >= 0 ? newDToxin : 0;

		int newEnergy = element.cell.energy - params.digestibleToxinDigestionCost;
		element.cell.energy = newEnergy >= 0 ? newEnergy : 0;
	}
}

__device__
void create_d_toxin(GridElement &element, CellModelParams &params) {
	int newDToxin = element.cell.dToxin + params.digestibleToxinGenerationRate;
	int maxDToxin = (1 << Cell::nbits_d_toxin) - 1;
	element.cell.dToxin = newDToxin <= maxDToxin ? newDToxin : maxDToxin;
}

__device__
int get_idx(int x, int y, int z, CellModelParams &params) {
	bool inBounds = x >= 0 && x < params.w
				 && y >= 0 && y < params.h
				 && z >= 0 && z < params.d;
	return inBounds ? x + y * params.w + z * params.w * params.h : -1;
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
	return pow(depth, 2.5);
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