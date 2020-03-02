#include <iostream>
#include <algorithm>
#include <curand_kernel.h>

#include "cell_model.cuh"
#include "param/simulation_parameters.h"
#include "operation/util.cuh"
#include "operation/movement.cuh"
#include "operation/growth.cuh"
#include "util/cuda_util.h"

// TODO: Refactor cuda functions to dedicated files:

// Cuda simulation functions:
__global__ void initialise_grid(GridElement *grid, ModelParameters params);
__global__ void update_cells(GridElement *grid, ModelParameters params);
__global__ void update_environment(GridElement *grid, ModelParameters params);
__global__ void update_interactions(GridElement *grid, ModelParameters params);
__device__ void iterate_cell(GridElement *grid, ModelParameters params, int idx, int x, int y, int z);

// GridElement initialisation functions:
__device__ void initialise_cell(Cell &cell, int idx, ModelParameters params);
__device__ void initialise_environment(GridElement &element, int idx, ModelParameters params);

// Cell operations:
__device__ void check_death(GridElement &element, ModelParameters &params);
__device__ void use_energy(GridElement &element, int y, ModelParameters &params);
__device__ void use_chem(GridElement &element, int y, ModelParameters &params);
__device__ void acquire_energy(GridElement &element, int y, ModelParameters &params);
__device__ void acquire_chem(GridElement &element, int x, int y, int z, ModelParameters &params);
__device__ void create_d_toxin(GridElement &element, ModelParameters &params);
__device__ void digest_d_toxin(GridElement &element, ModelParameters &params);

// Environment distributions:
__device__ double light_distribution(int y, int height);
__device__ double co2_distribution(int y, int height);
__device__ double chemical_distribution(int y, int height);

CellModel::CellModel(SimulationParameters params) :
		params(params),
		numBlocks(params.cuda.numBlocks),
		blockSize(params.cuda.blockSize) {
	int gridSize = params.model.w * params.model.h * params.model.d;
	std::cout << "BLOCKS: " << params.cuda.numBlocks << std::endl;
	cudaMallocManaged(&grid, gridSize * sizeof(GridElement));
	checkCudaError(cudaPeekAtLastError());
	initialise();
}

GridElement* CellModel::getHostGrid() {
	synchronizeData();
	return grid;
}

GridElement* CellModel::getDeviceGrid() {
	return grid;
}

SimulationParameters CellModel::getParams() {
	return params;
}

void CellModel::initialise() {
	initialise_grid<<<numBlocks, blockSize>>>(
		grid,
		params.model
	);
	checkCudaError(cudaPeekAtLastError());
}

void CellModel::simulate(int nIterations) {
	for (int i = 0; i < nIterations; i++) {
		// Cells:
		update_cells<<<numBlocks, blockSize>>>(grid, params.model);
		checkCudaError(cudaPeekAtLastError());

		// Environment
		update_environment<<<numBlocks, blockSize>>>(grid, params.model);
		checkCudaError(cudaPeekAtLastError());

		// Simulation:
		update_interactions<<<numBlocks, blockSize>>>(grid, params.model);
		checkCudaError(cudaPeekAtLastError());
	}
}

void CellModel::synchronizeData() {
	checkCudaError(cudaPeekAtLastError());
	cudaDeviceSynchronize();
	checkCudaError(cudaPeekAtLastError());
}

__global__
void initialise_grid(GridElement *grid, ModelParameters params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		grid[idx] = {};

		GridElement &element = grid[idx];
		curand_init(params.gridRandomSeed, idx, 0, &element.randState);
		curand_init(params.cellRandomSeed, idx, 0, &element.cell.randState);
		element.position.idx = idx;
		element.position.x = idx % params.w;
		element.position.y = (idx / params.w) % params.h;
		element.position.z = idx / (params.w * params.h);

		if (curand_uniform(&grid[idx].randState) < params.initialCellDensity) {
			initialise_cell(element.cell, idx, params);
		}

		initialise_environment(element, idx, params);
	}
}

__device__
void initialise_cell(Cell &cell, int idx, ModelParameters params) {
	cell.alive = true;
	cell.energy = 230;
	cell.chem = 230;
	cell.dToxin = 10;
	cell.ndToxin = 10;
	cell.parent_idx = idx;
}

__device__
void initialise_environment(GridElement &element, int idx, ModelParameters params) {
	curand_init(params.environmentRandomSeed, idx, 0, &element.environment.randState);
	if (curand_uniform(&element.environment.randState) < params.initialChemDensity) {
		int y = (idx / params.w) % params.h;
		element.environment.chem = chemical_distribution(y, params.h) * params.initialChemMax;
	}
	if (curand_uniform(&element.environment.randState) < params.initialNdToxinDensity) {
		element.environment.ndToxin = curand_uniform(&element.environment.randState) * params.initialNdToxinMax;
	}
	element.environment.dToxin = 0;
}

__global__
void update_cells(GridElement *grid, ModelParameters params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		int x = idx % params.w;
		int y = (idx / params.w) % params.h;
		int z = idx / (params.w * params.h);
		iterate_cell(grid, params, idx, x, y, z);
	}
}

__global__
void update_environment(GridElement *grid, ModelParameters params) {

}

__global__
void update_interactions(GridElement *grid, ModelParameters params) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		if (grid[idx].cell.alive) {
			// movement::execute(grid, grid[idx], params);
		} else {
			growth::execute(grid, grid[idx], params);
		}
	}
}

__device__
void iterate_cell(GridElement *grid, ModelParameters params, int idx, int x, int y, int z) {
	GridElement &element = grid[idx];

	if (element.cell.alive) {
		// check_death(element, params);
		use_energy(element, y, params);
		use_chem(element, y, params);
		acquire_energy(element, y, params);
		acquire_chem(element, x, y, z, params);
		create_d_toxin(element, params);
		digest_d_toxin(element, params);

		// Movement:
		// movement::prepare(grid, element, params);
		growth::prepare(grid, element, params);
	}
}

__device__
void check_death(GridElement &element, ModelParameters &params) {
	if (element.cell.has_subcell || element.cell.is_subcell) {
		return;
	}
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
void use_energy(GridElement &element, int y, ModelParameters &params) {
	int newEnergy = element.cell.energy - params.energyUsageRate;
	element.cell.energy = newEnergy > 0 ? newEnergy : 0;
}

__device__
void use_chem(GridElement &element, int y, ModelParameters &params) {
	int newChem = element.cell.chem - params.chemUsageRate;
	element.cell.chem = newChem > 0 ? newChem : 0;
}

__device__
void acquire_energy(GridElement &element, int y, ModelParameters &params) {
	int newEnergy = element.cell.energy;

	newEnergy += params.lightEnergyConversionRate * light_distribution(y, params.h);
	newEnergy += params.co2EnergyConversionRate * co2_distribution(y, params.h);

	int maxEnergy = (1 << Cell::nbits_energy) - 1;
	element.cell.energy = newEnergy <= maxEnergy ? newEnergy : maxEnergy;
}

__device__
void acquire_chem(GridElement &element, int x, int y, int z, ModelParameters &params) {
	int quantity = min(params.chemAcquisitionRate, element.environment.chem);
	
	int maxChem = (1 << Cell::nbits_chem) - 1;
	quantity = min(quantity, maxChem - element.cell.chem);

	element.cell.chem += quantity;
	element.environment.chem -= quantity;
}

__device__
void digest_d_toxin(GridElement &element, ModelParameters &params) {
	if (element.cell.dToxin >= params.dToxinDigestionThreshold) {
		int newDToxin = element.cell.dToxin - params.digestibleToxinDigestionRate;
		element.cell.dToxin = newDToxin >= 0 ? newDToxin : 0;

		int newEnergy = element.cell.energy - params.digestibleToxinDigestionCost;
		element.cell.energy = newEnergy >= 0 ? newEnergy : 0;
	}
}

__device__
void create_d_toxin(GridElement &element, ModelParameters &params) {
	int newDToxin = element.cell.dToxin + params.digestibleToxinGenerationRate;
	int maxDToxin = (1 << Cell::nbits_d_toxin) - 1;
	element.cell.dToxin = newDToxin <= maxDToxin ? newDToxin : maxDToxin;
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