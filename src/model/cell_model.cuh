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

	std::chrono::duration<double> durationCells;
	std::chrono::duration<double> durationBigCells;
	std::chrono::duration<double> durationInteractions;
	std::chrono::duration<double> durationPrepareGrowth;
	std::chrono::duration<double> durationGrowthInteractions;
	std::chrono::duration<double> durationEnvironment;

private:
	void initialise();
	void iterateRandomMovement();

	GridElement *grid;
	SimulationParameters params;
	dim3 blockSize;
	dim3 numBlocks;
	int iterations = 0;


};
