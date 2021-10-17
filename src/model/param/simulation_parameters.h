#pragma once

#include <string>

using std::string;

class Intervention;

struct CudaParameters {
	int blockSize = 1024;
	int numBlocks = 10;

	// TODO: Reduction buffer params
	// TODO: Seed parameter
};

struct OutputParameters {
	struct Video {
		bool enabled;
		int interval = 1000;
		string energyPath;
		string chemPath;
		string toxinPath;
		string genomePath;
		bool active = false;
	} video;

	struct Statistics {
		bool enabled;
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

	double lightIntensity;
	double temperature;
	double optimalTemperature;
	double functionalTemperatureRange;

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

	// Genomes:
	int genomeNum;
	int startingEnergy[10];
	int startingChem[10];
	int energyUsageRate[10];
	int chemUsageRate[10];

};

class SimulationParameters {
public:
	SimulationParameters();
	SimulationParameters(int w, int h, int d);

	CudaParameters cuda;
	OutputParameters output;
	ModelParameters model;

	int nInterventions = 0;
	Intervention **interventions;

	static SimulationParameters fromJson(string path);
	static SimulationParameters getDefaultParams();

	int gridSize();
};
