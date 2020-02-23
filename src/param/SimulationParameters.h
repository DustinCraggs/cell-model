#pragma once
#include <string>

using std::string;

struct CudaParams {
	int blockSize = 1024;
	int numBlocks = 10;

	// TODO: Reduction buffer params
	// TODO: Seed parameter
};

struct OutputParams {
	struct Video {
		int interval;
		string energyPath;
		string chemPath;
		string toxinPath;
	} video;

	struct Statistics {
		int interval;
		bool console;
		string csvPath;
	} statistics;
};

struct CellModelParams {
	// Initialisation params:
	int w, h, d;
	int iterations;
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
	int energyUsageRate;
	int chemUsageRate;
	int chemAcquisitionRate;
	int lightEnergyConversionRate;
	int co2EnergyConversionRate;
	int digestibleToxinGenerationRate;
	int digestibleToxinDigestionRate;
	int digestibleToxinDigestionCost;

	// Movement:
	float movementProbability;

	CudaParams cudaParams;

	CellModelParams(int w, int h, int d);

	static CellModelParams fromJson(string path, OutputParams *outputParams);
	static CellModelParams getDefaultParams(int w, int h, int d);
};

