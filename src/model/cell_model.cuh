#pragma once

#include <chrono>
#include "grid_element.h"
#include "param/simulation_parameters.h"

class CellModel {
public:
	CellModel(SimulationParameters params);
	void printCells();
	void simulate(int nIterations);
	void synchronizeData();

	GridElement* getHostGrid();
	GridElement* getDeviceGrid();
	SimulationParameters getParams();
	void updateParams(SimulationParameters params);
	void redistributeChemicals(double newChemDensity, bool invertDistribution);

	int runtimeCells;
	int runtimeBigCells;
	int runtimeInteractions;
	int runtimePrepareGrowth;
	int runtimeGrowthInteractions;
	int runtimeEnvironment;

	int totalCellsRuntime = 0;
	int totalBigCellsRuntime = 0;	
	int totalInteractionsRuntime = 0;
	int totalPrepareGrowthRuntime = 0;
	int totalGrowthInteractionsRuntime = 0;
	int totalEnvironmentRuntime = 0;

private:
	void initialise();
	void iterateRandomMovement();

	GridElement *grid;
	SimulationParameters params;
	dim3 blockSize;
	dim3 numBlocks;
	int iterations = 0;


};
