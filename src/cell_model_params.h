#pragma once

struct CudaParams {
	dim3 blockSize = 1024;
	dim3 numBlocks = 10;

	// TODO: Reduction buffer params
	// TODO: Seed parameter
};

struct CellModelParams {
	// Initialisation params:
	int w, h, d;
	int gridSize;
	float initialCellDensity;
	float initialChemDensity;
	float initialChemMax;
	float initialNdToxinDensity;
	float initialNdToxinMax;

	// Death:
	int energySurvivalThreshold;
	int chemSurvivalThreshold;
	int dToxinDeathThreshold;
	int dToxinDigestionThreshold;
	int ndToxinDeathThreshold;

	// Resource rates:
	float energyUsageRate;
	float chemUsageRate;
	float lightEnergyConversionRate;
	float co2EnergyConversionRate;
	float digestibleToxinGenerationRate;
	float digestibleToxinDigestionRate;
	float digestibleToxinDigestionCost;

	int chemAcquisitionRate;

	// Movement:
	float movementProbability;

	CudaParams cudaParams;

	CellModelParams(int w, int h, int d) {
		this->w = w;
		this->h = h;
		this->d = d;
		this->gridSize = w * h * d;
	}
};