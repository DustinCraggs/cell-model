#include <iostream>
#include <stdio.h>

#include "cell_model.cuh"

// TODO: Const things

int main() {
	CellModelParams params(5, 5);
	params.initialDensity = 0.5;
	params.survivalThreshold = 3;
	params.energyLossP = 0.9;
	params.gatherLightEnergyP = 0.5;

	// Cuda parameters:
	params.cudaParams.blockSize = 1024;
	params.cudaParams.numBlocks = 10;
	
	CellModel model(params);
	model.printCells();

	for (int i = 0; i < 10; i++) {
		model.simulate(1);
		model.printCells();
		std::cout << std::endl << std::endl;
	}
}