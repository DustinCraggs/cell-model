#include <iostream>
#include <chrono>
// #include <stdio.h>

#include "cell_model.cuh"

// TODO: Const things
#define CURRENT_TIME std::chrono::steady_clock::now()
#define TIME_DIFF(T0, T1) std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count()

void printStatistics(CellModel model) {
	auto t0 = CURRENT_TIME;
	int nLiving = model.numberOfLivingCells();
	std::cout << TIME_DIFF(t0, CURRENT_TIME) << std::endl;

	std::cout << "N living cells: " << nLiving << std::endl;
	if (nLiving != 0) {
		std::cout << "Enrg avg.: " << model.totalEnergy()/nLiving << std::endl;
	}
	std::cout << std::endl;
}

int main() {
	CellModelParams params(100, 100, 10);
	std::cout << params.nCells << std::endl;
	params.initialDensity = 1.1;
	params.survivalThreshold = 3;
	params.energyLossProb = 0.35;
	params.gatherLightEnergyProb = 0.3;

	// Cuda parameters:
	params.cudaParams.blockSize = 1024;
	params.cudaParams.numBlocks = 10;
	
	CellModel model(params);
	// model.synchronizeData();
	// model.printCells();

	// printStatistics(model);
	for (int i = 0; i < 30; i++) {
		model.simulate(100);
		printStatistics(model);

		// model.synchronizeData();
		// model.printCells();
	}
}