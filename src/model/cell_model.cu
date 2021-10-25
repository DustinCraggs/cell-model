#include <iostream>
#include <algorithm>
#include <chrono>
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
__global__ void update_cells(GridElement *grid, ModelParameters params, int iterations);
__global__ void update_big_cells(GridElement *grid, ModelParameters params, int iterations);
__global__ void prepare_growth(GridElement *grid, ModelParameters params, int iterations);
__global__ void update_environment(GridElement *grid, ModelParameters &params);
__global__ void update_interactions(GridElement *grid, ModelParameters params, int iterations);
__global__ void update_growth_interactions(GridElement *grid, ModelParameters params, int iterations);
__device__ void iterate_cell(GridElement *grid, ModelParameters params, int idx, int x, int y, int z, int iteration);

// Interventions
__global__ void add_new_chemicals(GridElement *grid, ModelParameters params,
	double newChemDensity, bool invertDistribution);

// GridElement initialisation functions:
__device__ void initialise_cell(Cell &cell, int idx, ModelParameters params, float randNum);
__device__ void initialise_environment(GridElement &element, int idx, ModelParameters &params);
__device__ void add_chem_to_element(GridElement &element, ModelParameters &params, double density,
	int maxAddedChem, bool invertDistribution);

// Cell operations:
__device__ void check_death(GridElement &element, ModelParameters &params);
__device__ void use_energy(GridElement &element, int y, ModelParameters &params);
__device__ void use_chem(GridElement &element, int y, ModelParameters &params);
__device__ void acquire_energy(GridElement &element, int y, ModelParameters &params);
__device__ void acquire_chem(GridElement &element, int x, int y, int z, ModelParameters &params);
__device__ void create_d_toxin(GridElement &element, ModelParameters &params);
__device__ void digest_d_toxin(GridElement &element, ModelParameters &params);


// Environment distributions:
__device__ double energy_efficiency(double temperature, double optimalTemperature,
	double functionalTemperatureRange);
__device__ double light_distribution(int y, int height);
__device__ double co2_distribution(int y, int height);
__device__ double chemical_distribution(int y, int height, bool invertDistribution);

CellModel::CellModel(SimulationParameters params) :
		params(params),
		numBlocks(params.cuda.numBlocks),
		blockSize(params.cuda.blockSize) {
	int gridSize = params.model.w * params.model.h * params.model.d;
	std::cout << "BLOCKS: " << params.cuda.numBlocks << std::endl;
	cudaMallocManaged(&grid, gridSize * sizeof(GridElement));
	std::cout << gridSize << std::endl;
	std::cout << sizeof(GridElement) << std::endl;
	std::cout << gridSize * sizeof(GridElement) << std::endl;
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

void CellModel::updateParams(SimulationParameters params) {
	this->params = params;
}

void CellModel::redistributeChemicals(double newChemDensity, bool invertDistribution) {
	add_new_chemicals<<<numBlocks, blockSize>>>(grid, params.model, newChemDensity, invertDistribution);
	checkCudaError(cudaPeekAtLastError());
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
		// Cells
		auto start = std::chrono::high_resolution_clock::now();
		update_cells<<<numBlocks, blockSize>>>(grid, params.model, iterations);
		auto stop = std::chrono::high_resolution_clock::now();
		runtimeCells = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		checkCudaError(cudaPeekAtLastError());

		start = std::chrono::high_resolution_clock::now();
		update_big_cells<<<numBlocks, blockSize>>>(grid, params.model, iterations);
		stop = std::chrono::high_resolution_clock::now();
		runtimeBigCells = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		checkCudaError(cudaPeekAtLastError());

		start = std::chrono::high_resolution_clock::now();
		update_interactions<<<numBlocks, blockSize>>>(grid, params.model, iterations);
		stop = std::chrono::high_resolution_clock::now();
		runtimeInteractions = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		checkCudaError(cudaPeekAtLastError());

		// Growth:
		start = std::chrono::high_resolution_clock::now();
		prepare_growth<<<numBlocks, blockSize>>>(grid, params.model, iterations);
		stop = std::chrono::high_resolution_clock::now();
		runtimePrepareGrowth = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		checkCudaError(cudaPeekAtLastError());

		start = std::chrono::high_resolution_clock::now();
		update_growth_interactions<<<numBlocks, blockSize>>>(grid, params.model, iterations);
		stop = std::chrono::high_resolution_clock::now();
		runtimeGrowthInteractions = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		checkCudaError(cudaPeekAtLastError());

		// Environment
		start = std::chrono::high_resolution_clock::now();
		update_environment<<<numBlocks, blockSize>>>(grid, params.model);
		stop = std::chrono::high_resolution_clock::now();
		runtimeEnvironment = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		checkCudaError(cudaPeekAtLastError());

		// Simulation:
		this->iterations++;
	}

}

__global__
void update_big_cells(GridElement *grid, ModelParameters params, int iterations) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		// Large cells:
		if (grid[idx].cell.alive && grid[idx].cell.has_subcell && !grid[idx].cell.is_subcell) {
			growth::checkDeathAndDistributeResources(grid, grid[idx], params);
		}
		// growth::prepare(grid, grid[idx], params);
	}
}

__global__
void prepare_growth(GridElement *grid, ModelParameters params, int iterations) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		if (grid[idx].cell.alive) {
			growth::prepare(grid, grid[idx], params, iterations);
		}
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
		element.canGrow = false;

		float randNum = curand_uniform(&grid[idx].randState);

		int maxGenomeNum = 10;

		// Different modes for density (based on genome num or not)
		if(randNum < (params.initialCellDensity)*(float(params.genomeNum)/maxGenomeNum)) {
		// if(randNum < params.initialCellDensity) {
			initialise_cell(element.cell, idx, params, randNum);
		}

		initialise_environment(element, idx, params);
	}
}

__device__
void initialise_cell(Cell &cell, int idx, ModelParameters params, float randNum) {

	cell.alive = true;
	cell.energy = 230;
	cell.chem = 230;
	cell.dToxin = 10;
	cell.ndToxin = 10;
	cell.nextCellOffset = {0};
	cell.parent_idx = idx;
	cell.is_subcell = false;
	cell.has_subcell = false;

	// Keep as float since it needs to be for fraction calculation.
	float genomeNum = params.genomeNum;

	int maxGenomeNum = 10;

	for(int i = 1; i <= genomeNum; i++) {

		// Different modes for density (based on genome num or not)
		if(randNum <= (params.initialCellDensity*((1/genomeNum)*i))*(float(genomeNum/maxGenomeNum))) {
		// if(randNum <= params.initialCellDensity*((1/genomeNum)*i)) {

			cell.genome = i;
			break;

		}

	}

	int genomeIndex = cell.genome - 1;

	cell.energy = params.startingEnergy[genomeIndex];
	cell.chem = params.startingChem[genomeIndex];
	cell.energyUsageRate = params.energyUsageRate[genomeIndex];
	cell.chemUsageRate = params.chemUsageRate[genomeIndex];

}

__device__
void initialise_environment(GridElement &element, int idx, ModelParameters &params) {
	curand_init(params.environmentRandomSeed, idx, 0, &element.environment.randState);

	add_chem_to_element(element, params, params.initialChemDensity, params.initialChemMax, false);

	if (curand_uniform(&element.environment.randState) < params.initialNdToxinDensity) {
		element.environment.ndToxin = curand_uniform(&element.environment.randState) * params.initialNdToxinMax;
	}
	element.environment.dToxin = 0;
}

__device__
void add_chem_to_element(GridElement &element, ModelParameters &params,
		double density, int maxAddedChem, bool invertDistribution) {
	if (curand_uniform(&element.environment.randState) < density) {
		int y = (element.position.idx / params.w) % params.h;
		int quantity = chemical_distribution(y, params.h, invertDistribution) * maxAddedChem;

		int maxChem = (1 << Environment::nbits_chem) - 1;
		element.environment.chem = min(maxChem, element.environment.chem + quantity);
	}
}

__global__
void update_cells(GridElement *grid, ModelParameters params, int iterations) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		int x = idx % params.w;
		int y = (idx / params.w) % params.h;
		int z = idx / (params.w * params.h);
		iterate_cell(grid, params, idx, x, y, z, iterations);
	}
}

__global__
void update_environment(GridElement *grid, ModelParameters &params) {

}

__global__
void update_growth_interactions(GridElement *grid, ModelParameters params, int iterations) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		if (grid[idx].cell.alive) {
			// movement::execute(grid, grid[idx], params);
		} else {
			// Large cells:
		}
		growth::execute(grid, grid[idx], params);
		// if (grid[idx].cell.has_subcell && !grid[idx].cell.is_subcell) {
		// 	growth::checkDeathAndDistributeResources(grid, grid[idx], params);
		// }
	}
}

__global__
void update_interactions(GridElement *grid, ModelParameters params, int iterations) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	
	// Stride loop:
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		if (grid[idx].cell.alive) {
			movement::execute(grid, grid[idx], params);
		} else {
			// growth::execute(grid, grid[idx], params);
		}
	}
}

// INTERVENTIONS:
__global__
void add_new_chemicals(GridElement *grid, ModelParameters params, double newChemDensity,
		bool invertDistribution) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int gridSize = params.w * params.h * params.d;
	for (int idx = tid; idx < gridSize; idx += blockDim.x * gridDim.x) {
		GridElement &element = grid[idx];
		add_chem_to_element(element, params, newChemDensity, params.initialChemMax, invertDistribution);
	}
}

__device__
void iterate_cell(GridElement *grid, ModelParameters params, int idx, int x, int y, int z, int iteration) {
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
		movement::prepare(grid, element, params);
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
	}
}

__device__
void use_energy(GridElement &element, int y, ModelParameters &params) {
	int newEnergy = element.cell.energy - element.cell.energyUsageRate;
	element.cell.energy = newEnergy > 0 ? newEnergy : 0;
}

__device__
void use_chem(GridElement &element, int y, ModelParameters &params) {
	int newChem = element.cell.chem - element.cell.chemUsageRate;
	element.cell.chem = newChem > 0 ? newChem : 0;
}

__device__
void acquire_energy(GridElement &element, int y, ModelParameters &params) {
	int newEnergy = 0;

	newEnergy += params.lightEnergyConversionRate
		* light_distribution(y, params.h)
		* params.lightIntensity;
	newEnergy += params.co2EnergyConversionRate * co2_distribution(y, params.h);
	newEnergy *= energy_efficiency(params.temperature,
		params.optimalTemperature, params.functionalTemperatureRange);

	newEnergy += element.cell.energy;

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
double energy_efficiency(double temperature, double optimalTemperature,
		double functionalTemperatureRange) {
	double variation = abs(temperature - optimalTemperature);
	double variationRatio = 1 - (variation / functionalTemperatureRange);
	return variationRatio > 0.0 ? variationRatio : 0.0;
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
double chemical_distribution(int y, int height, bool invertDistribution) {
	float depth = y/(float)(height - 1);
	if (invertDistribution) {
		depth = 1 - depth;
	}
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