#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

#include "cell_model_params.h"

using nlohmann::json;

CellModelParams::CellModelParams(int w, int h, int d) {
	this->w = w;
	this->h = h;
	this->d = d;
	this->gridSize = w * h * d;
}

CellModelParams CellModelParams::fromJson(std::string path, OutputParams *outputParams) {
	json jsonParams = json::parse(std::ifstream(path));

	if (!jsonParams.contains("model")) {
		// No model parameters in params
		return getDefaultParams(10, 10, 10);
	}

	json model = jsonParams["model"];
	int w = model.value("width", 10);
	int h = model.value("height", 10);
	int d = model.value("depth", 10);
	CellModelParams params = getDefaultParams(w, h, d);
	params.iterations = model.value("iterations", params.iterations);

	// Initialisation:
	params.initialCellDensity = model.value("initialCellDensity", params.initialCellDensity);
	params.initialChemDensity = model.value("initialChemDensity", params.initialChemDensity);
	params.initialChemMax = model.value("initialChemMax", params.initialChemMax);
	params.initialNdToxinDensity = model.value("initialNdToxinDensity", params.initialNdToxinDensity);
	params.initialNdToxinMax = model.value("initialNdToxinMax", params.initialNdToxinMax);

	// Thresholds:
	params.energySurvivalThreshold = model.value("energySurvivalThreshold", params.energySurvivalThreshold);
	params.chemSurvivalThreshold = model.value("chemSurvivalThreshold", params.chemSurvivalThreshold);
	params.dToxinDeathThreshold = model.value("dToxinDeathThreshold", params.dToxinDeathThreshold);
	params.dToxinDigestionThreshold = model.value("dToxinDigestionThreshold", params.dToxinDigestionThreshold);
	params.ndToxinDeathThreshold = model.value("ndToxinDeathThreshold", params.ndToxinDeathThreshold);

	// Energy rates:
	params.energyUsageRate = model.value("energyUsageRate", params.energyUsageRate);
	params.chemUsageRate = model.value("chemUsageRate", params.chemUsageRate);
	params.chemAcquisitionRate = model.value("chemAcquisitionRate", params.chemAcquisitionRate);
	params.lightEnergyConversionRate = model.value("lightEnergyConversionRate", params.lightEnergyConversionRate);
	params.co2EnergyConversionRate = model.value("co2EnergyConversionRate", params.co2EnergyConversionRate);
	params.digestibleToxinGenerationRate = model.value("digestibleToxinGenerationRate", params.digestibleToxinGenerationRate);
	params.digestibleToxinDigestionRate = model.value("digestibleToxinDigestionRate", params.digestibleToxinDigestionRate);
	params.digestibleToxinDigestionCost = model.value("digestibleToxinDigestionCost", params.digestibleToxinDigestionCost);

	// Movement:
	params.movementProbability = model.value("movementProbability", params.movementProbability);

	// Cuda parameters:
	if (jsonParams.contains("cuda")) {
		json cudaJsonParams = jsonParams["cuda"];
		params.cudaParams.blockSize = cudaJsonParams.value("blockSize", params.cudaParams.blockSize);
		params.cudaParams.numBlocks = cudaJsonParams.value("numBlocks", params.cudaParams.numBlocks);
	}

	if (jsonParams.contains("output")) {
		json outputJsonParams = jsonParams["output"];
		if (outputJsonParams.contains("video")) {
			json videoJsonParams = outputJsonParams["video"];
			(*outputParams).video.interval = videoJsonParams.value("interval", 1);
			(*outputParams).video.energyPath = videoJsonParams.value("energy", "");
			(*outputParams).video.chemPath = videoJsonParams.value("chemical", "");
			(*outputParams).video.toxinPath = videoJsonParams.value("toxin", "");
		}
	}

	return params;
}

CellModelParams CellModelParams::getDefaultParams(int w, int h, int d) {
	auto params = CellModelParams(w, h, d);
	params.iterations = 100;

	// Initialisation:
	params.initialCellDensity = 0.1;
	params.initialChemDensity = 0.3;
	params.initialChemMax = 255;
	params.initialNdToxinDensity = 0.05;
	params.initialNdToxinMax = 150;

	// Thresholds:
	params.energySurvivalThreshold = 20;
	params.chemSurvivalThreshold = 1;
	params.dToxinDeathThreshold = 254;
	params.dToxinDigestionThreshold = 100;
	params.ndToxinDeathThreshold = 254;

	// Energy rates:
	params.energyUsageRate = 15;
	params.chemUsageRate = 1;
	params.chemAcquisitionRate = 30;
	params.lightEnergyConversionRate = 10;
	params.co2EnergyConversionRate = 15;
	params.digestibleToxinGenerationRate = 1;
	params.digestibleToxinDigestionRate = 1;
	params.digestibleToxinDigestionCost = 1;

	// Movement:
	params.movementProbability = 0.1;

	// Cuda parameters:
	params.cudaParams.blockSize = 1024;
	params.cudaParams.numBlocks = 10;

	return params;
}