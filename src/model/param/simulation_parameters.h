#pragma once

#include <string>

using std::string;

struct CudaParameters {
	int blockSize = 1024;
	int numBlocks = 10;

	// TODO: Reduction buffer params
	// TODO: Seed parameter
};

struct OutputParameters {
	struct Video {
		int interval;
		string energyPath;
		string chemPath;
		string toxinPath;
	} video;

	struct Statistics {
		int interval;
		bool console;
		string file;
	} statistics;

	struct Runtime {
		string file;
	} runtime;
};

struct ModelParameters {
	// Initialisation params:
	int w, h, d;
	int iterations;

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

	// Growth:
	int maxCellSize;
	int growthCost;
	int growthThreshold;

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

	// RNG:
	int cellRandomSeed;
	int environmentRandomSeed;
	int gridRandomSeed;
};

class SimulationParameters {
public:
	SimulationParameters();
	SimulationParameters(int w, int h, int d);

	CudaParameters cuda;
	OutputParameters output;
	ModelParameters model;

	static SimulationParameters fromJson(string path);
	static SimulationParameters getDefaultParams();

	int gridSize();

};
