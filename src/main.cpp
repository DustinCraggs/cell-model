#include <iostream>
#include <chrono>

#include "model/cell_model_driver.h"

// TODO: Const things

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cerr << "Error: no path to configuration file provided" << std::endl;
		std::cerr << "E.g: simulate experiment/config/example.json" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	CellModelDriver driver(argv[1]);
	driver.run();
}



// #define CURRENT_TIME std::chrono::steady_clock::now()
// #define TIME_DIFF(T0, T1) std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count()

// void printStatisticsCsvRow(CellModel model, int iterationNumber) {
// 	int nLiving = model.numberOfLivingCells();
// 	float averageCellEnergy = nLiving > 0 ? model.totalCellEnergy()/nLiving : 0;
// 	float averageCellChem = nLiving > 0 ? model.totalCellChem()/nLiving : 0;
// 	float averageCellToxin = nLiving > 0 ? model.totalCellToxin()/nLiving : 0;

// 	// Environment:
// 	float totalEnvChem = model.totalEnvChem();
// 	float totalEnvToxin = model.totalEnvToxin();

// 	printf("%d,%d,%f,%f,%f,%f,%f\n", iterationNumber, nLiving, averageCellEnergy, averageCellChem,
// 		averageCellToxin, totalEnvChem, totalEnvToxin);
// }

// void printStatistics(CellModel model) {
// 	int nLiving = model.numberOfLivingCells();

// 	std::cerr << "N living cells: " << nLiving << std::endl;
// 	if (nLiving != 0) {
// 		std::cerr << "Enrg avg.: " << model.totalCellEnergy()/nLiving << std::endl;
// 	}
// 	std::cerr << std::endl;
// }

// int main(int argc, char **argv) {
// 	// // CellModelParameters params(8, 8, 1);
// 	// CellModelParameters params(128, 128, 1);
// 	// // CellModelParameters params(256, 256, 1);
// 	// // Initialisation:
// 	// params.initialCellDensity = 0.1;
// 	// params.initialChemDensity = 0.3;
// 	// params.initialChemMax = 255;
// 	// params.initialNdToxinDensity = 0.05;
// 	// params.initialNdToxinMax = 150;

// 	// // Thresholds:
// 	// params.energySurvivalThreshold = 20;
// 	// params.chemSurvivalThreshold = 1;
// 	// params.dToxinDeathThreshold = 254;
// 	// params.dToxinDigestionThreshold = 100;
// 	// params.ndToxinDeathThreshold = 254;

// 	// // Energy rates:
// 	// params.energyUsageRate = 15;
// 	// params.chemUsageRate = 1;
// 	// params.lightEnergyConversionRate = 10;
// 	// params.co2EnergyConversionRate = 15;
// 	// params.digestibleToxinGenerationRate = 1;
// 	// params.digestibleToxinDigestionRate = 1;
// 	// params.digestibleToxinDigestionCost = 1;

// 	// params.chemAcquisitionRate = 30;

// 	// // Movement:
// 	// params.movementProbability = 0.1;

// 	// // Cuda parameters:
// 	// params.cudaParams.blockSize = 1024;
// 	// params.cudaParams.numBlocks = 10;
	
// 	auto params = CellModelParameters::fromJson(argv[1]);
// 	CellModel model(params);
// 	model.synchronizeData();

// 	// model.writeFrame(0);
// 	std::cerr << "Iteration: " << 0 << std::endl;
// 	printStatistics(model);
// 	// printStatisticsCsvRow(model, 0);

// 	for (int i = 0; i < 500; i++) {
// 		model.simulate(1);
// 		model.synchronizeData();
// 		if ((i + 1) % 50 == 0) {
// 			std::cerr << "Iteration: " << i + 1 << std::endl;
// 			printStatistics(model);
// 		}
// 		if ((i + 1) % 10 == 0) {
// 			printStatisticsCsvRow(model, i);
// 		}
// 		// model.writeFrame(0);
// 	}
// }

// auto t0 = CURRENT_TIME;
// std::cerr << TIME_DIFF(t0, CURRENT_TIME) << std::endl;