#include <iostream>
#include <algorithm>
#include <curand_kernel.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>

#include "cell_model.cuh"
#include "cell_model_params.h"
#include "cuda_util.h"
#include "cell/movement.cuh"

#include <curand_kernel.h>

#define VIDEO_TYPE_TOXIN

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

void CellModel::writeFrame(int dIdx) {
	int nPixels = params.w * params.h;
	unsigned char* framebuffer = new unsigned char[nPixels * 3];

	for (int i = 0; i < nPixels; i++) {
		GridElement element = grid[i + dIdx*nPixels];
		#ifdef VIDEO_TYPE_ENERGY
			if (element.cell.alive) {
				framebuffer[i*3] = 0.5 * (255 - element.cell.energy);
				framebuffer[i*3 + 1] = 0;
				framebuffer[i*3 + 2] = 0.8 * (element.cell.energy);
			} else {
				framebuffer[i*3] = 0;
				framebuffer[i*3 + 1] = element.environment.chem * 0.6;
				framebuffer[i*3 + 2] = 0;
			}
		#endif
		#ifdef VIDEO_TYPE_CHEM
			if (element.cell.alive) {
				framebuffer[i*3] = 0.5 * (255 - element.cell.chem);
				framebuffer[i*3 + 1] = 0;
				framebuffer[i*3 + 2] = 0.8 * (element.cell.chem);
			} else {
				framebuffer[i*3] = 0;
				framebuffer[i*3 + 1] = element.environment.chem * 0.6;
				framebuffer[i*3 + 2] = 0;
			}
		#endif
		#ifdef VIDEO_TYPE_TOXIN
			if (element.cell.alive) {
				framebuffer[i*3] = 0.4 * (255 - element.cell.dToxin);
				framebuffer[i*3 + 1] = 0.4 * (element.cell.dToxin);
				framebuffer[i*3 + 2] = 0.4 * (element.cell.dToxin);
			} else {
				framebuffer[i*3] = (element.environment.ndToxin + element.environment.dToxin) * 0.2;
				framebuffer[i*3 + 1] = (element.environment.ndToxin + element.environment.dToxin) * 0.2;
				framebuffer[i*3 + 2] = 0;
			}
		#endif
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

struct CellChem {
	__device__
	double operator()(const GridElement& g) const {
		return g.cell.alive ? g.cell.chem : 0.0;
	}
};

struct CellToxin {
	__device__
	double operator()(const GridElement& g) const {
		return g.cell.alive ? g.cell.dToxin + g.cell.ndToxin : 0.0;
	}
};

struct EnvChem {
	__device__
	double operator()(const GridElement& g) const {
		return g.environment.chem;
	}
};

struct EnvToxin {
	__device__
	double operator()(const GridElement& g) const {
		return g.environment.dToxin + g.environment.ndToxin;
	}
};

int CellModel::numberOfLivingCells() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellOccupied(), 0, thrust::plus<int>());
}

double CellModel::totalCellEnergy() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellEnergy(), 0.0, thrust::plus<double>());
}

double CellModel::totalCellChem() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellChem(), 0.0, thrust::plus<double>());
}

double CellModel::totalCellToxin() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellToxin(), 0.0, thrust::plus<double>());
}

double CellModel::totalEnvChem() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, EnvChem(), 0.0, thrust::plus<double>());
}

double CellModel::totalEnvToxin() {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, EnvToxin(), 0.0, thrust::plus<double>());
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