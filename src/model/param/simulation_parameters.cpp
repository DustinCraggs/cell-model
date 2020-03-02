#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "simulation_parameters.h"

using nlohmann::json;

void cudaParametersFromJson(CudaParameters &cuda, json jsonParams);
void outputParametersFromJson(OutputParameters &output, json jsonParams);
void modelParametersFromJson(ModelParameters &model, json jsonParams);

SimulationParameters::SimulationParameters() {

}

SimulationParameters::SimulationParameters(int w, int h, int d) {
	this->model.w = w;
	this->model.h = h;
	this->model.d = d;
}

int SimulationParameters::gridSize() {
	return model.w * model.h * model.d;
}

SimulationParameters SimulationParameters::fromJson(std::string path) {
	json jsonParams = json::parse(std::ifstream(path));
	SimulationParameters simulationParameters = getDefaultParams();

	if (jsonParams.contains("cuda")) {
		cudaParametersFromJson(simulationParameters.cuda, jsonParams["cuda"]);
	}

	if (jsonParams.contains("output")) {
		outputParametersFromJson(simulationParameters.output, jsonParams["output"]);
	}
	
	if (jsonParams.contains("model")) {
		modelParametersFromJson(simulationParameters.model, jsonParams["model"]);
	}

	return simulationParameters;
}

void cudaParametersFromJson(CudaParameters &cuda, json jsonParams) {
	cuda.blockSize = jsonParams.value("blockSize", cuda.blockSize);
	cuda.numBlocks = jsonParams.value("numBlocks", cuda.numBlocks);
}

void outputParametersFromJson(OutputParameters &output, json jsonParams) {
	if (jsonParams.contains("video")) {
		json videoJsonParams = jsonParams["video"];
		output.video.interval = videoJsonParams.value("interval", 1);
		output.video.energyPath = videoJsonParams.value("energy", "");
		output.video.chemPath = videoJsonParams.value("chemical", "");
		output.video.toxinPath = videoJsonParams.value("toxin", "");
	}
	if (jsonParams.contains("statistics")) {
		json statisticsJsonParams = jsonParams["statistics"];
		output.statistics.interval = statisticsJsonParams.value("interval", 1);
		output.statistics.file = statisticsJsonParams.value("file", "");
		output.statistics.console = statisticsJsonParams.value("console", false);
	}
	if (jsonParams.contains("runtime")) {
		json runtimeJsonParams = jsonParams["runtime"];
		output.runtime.file = runtimeJsonParams.value("file", "");
	}
}

void modelParametersFromJson(ModelParameters &model, json jsonParams) {
	model.w = jsonParams.value("width", model.w);
	model.h = jsonParams.value("height", model.h);
	model.d = jsonParams.value("depth", model.d);
	model.iterations = jsonParams.value("iterations", model.iterations);
	// Initialisation:
	model.initialCellDensity = jsonParams.value("initialCellDensity", model.initialCellDensity);
	model.initialChemDensity = jsonParams.value("initialChemDensity", model.initialChemDensity);
	model.initialChemMax = jsonParams.value("initialChemMax", model.initialChemMax);
	model.initialNdToxinDensity = jsonParams.value("initialNdToxinDensity", model.initialNdToxinDensity);
	model.initialNdToxinMax = jsonParams.value("initialNdToxinMax", model.initialNdToxinMax);
	// Thresholds:
	model.energySurvivalThreshold = jsonParams.value("energySurvivalThreshold", model.energySurvivalThreshold);
	model.chemSurvivalThreshold = jsonParams.value("chemSurvivalThreshold", model.chemSurvivalThreshold);
	model.dToxinDeathThreshold = jsonParams.value("dToxinDeathThreshold", model.dToxinDeathThreshold);
	model.dToxinDigestionThreshold = jsonParams.value("dToxinDigestionThreshold", model.dToxinDigestionThreshold);
	model.ndToxinDeathThreshold = jsonParams.value("ndToxinDeathThreshold", model.ndToxinDeathThreshold);
	// Growth:
	model.maxCellSize = jsonParams.value("maxCellSize", model.maxCellSize);
	model.growthCost = jsonParams.value("growthCost", model.growthCost);
	model.growthThreshold = jsonParams.value("growthThreshold", model.growthThreshold);
	// Energy rates:
	model.energyUsageRate = jsonParams.value("energyUsageRate", model.energyUsageRate);
	model.chemUsageRate = jsonParams.value("chemUsageRate", model.chemUsageRate);
	model.chemAcquisitionRate = jsonParams.value("chemAcquisitionRate", model.chemAcquisitionRate);
	model.lightEnergyConversionRate = jsonParams.value("lightEnergyConversionRate", model.lightEnergyConversionRate);
	model.co2EnergyConversionRate = jsonParams.value("co2EnergyConversionRate", model.co2EnergyConversionRate);
	model.digestibleToxinGenerationRate = jsonParams.value("digestibleToxinGenerationRate", model.digestibleToxinGenerationRate);
	model.digestibleToxinDigestionRate = jsonParams.value("digestibleToxinDigestionRate", model.digestibleToxinDigestionRate);
	model.digestibleToxinDigestionCost = jsonParams.value("digestibleToxinDigestionCost", model.digestibleToxinDigestionCost);
	// Movement:
	model.movementProbability = jsonParams.value("movementProbability", model.movementProbability);
	//RNG:
	if (jsonParams.contains("randomSeed")) {
		json randomSeedJson = jsonParams["randomSeed"];
		model.cellRandomSeed = randomSeedJson.value("cell", model.cellRandomSeed);
		model.environmentRandomSeed = randomSeedJson.value("environment", model.cellRandomSeed);
		model.gridRandomSeed = randomSeedJson.value("grid", model.cellRandomSeed);
	}
}

SimulationParameters SimulationParameters::getDefaultParams() {
	auto params = SimulationParameters(100, 100, 50);

	// Model:
	auto &model = params.model;
	model.iterations = 100;
	// Initialisation:
	model.initialCellDensity = 0.1;
	model.initialChemDensity = 0.3;
	model.initialChemMax = 255;
	model.initialNdToxinDensity = 0.05;
	model.initialNdToxinMax = 150;
	// Thresholds:
	model.energySurvivalThreshold = 20;
	model.chemSurvivalThreshold = 1;
	model.dToxinDeathThreshold = 254;
	model.dToxinDigestionThreshold = 100;
	model.ndToxinDeathThreshold = 254;
	// Growth:
	model.maxCellSize = 10;
	model.growthCost = 20;
	model.growthThreshold = 250;
	// Energy rates:
	model.energyUsageRate = 15;
	model.chemUsageRate = 1;
	model.chemAcquisitionRate = 30;
	model.lightEnergyConversionRate = 10;
	model.co2EnergyConversionRate = 15;
	model.digestibleToxinGenerationRate = 1;
	model.digestibleToxinDigestionRate = 1;
	model.digestibleToxinDigestionCost = 1;
	// Movement:
	model.movementProbability = 0.1;
	// RNG:
	model.cellRandomSeed = 1000;
	model.environmentRandomSeed = 2000;
	model.gridRandomSeed = 3000;

	// Cuda parameters:
	auto &cuda = params.cuda;
	cuda.blockSize = 1024;
	cuda.numBlocks = 10;

	return params;
}