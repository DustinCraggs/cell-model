#include <iostream>
#include <chrono>

#include "cell_model.cuh"

// TODO: Const things

#define CURRENT_TIME std::chrono::steady_clock::now()
#define TIME_DIFF(T0, T1) std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count()

void printStatistics(CellModel model) {
	auto t0 = CURRENT_TIME;
	int nLiving = model.numberOfLivingCells();
	// std::cerr << TIME_DIFF(t0, CURRENT_TIME) << std::endl;

	std::cerr << "N living cells: " << nLiving << std::endl;
	if (nLiving != 0) {
		std::cerr << "Enrg avg.: " << model.totalEnergy()/nLiving << std::endl;
	}
	std::cerr << std::endl;
}

int main(int argc, char **argv) {
	CellModelParams params(128, 128, 1);
	params.initialDensity = 0.2;
	params.survivalThreshold = 20;
	params.energyUsageRate = 20;
	params.lightEnergyConversionRate = 15;
	params.co2EnergyConversionRate = 10;
	params.movementProbability = 0.01;

	// Cuda parameters:
	params.cudaParams.blockSize = 1024;
	params.cudaParams.numBlocks = 10;
	
	CellModel model(params);
	model.synchronizeData();

	model.writeFrame(0);
	std::cerr << "Iteration: " << 0 << std::endl;
	printStatistics(model);

	for (int i = 0; i < 500; i++) {
		model.simulate(1);
		if ((i + 1) % 50 == 0) {
			std::cerr << "Iteration: " << i + 1 << std::endl;
			printStatistics(model);
		}
		model.synchronizeData();
		model.writeFrame(0);
	}
}