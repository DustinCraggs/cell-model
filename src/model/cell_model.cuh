#pragma once

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
private:
	void initialise();
	void iterateRandomMovement();

	GridElement *grid;
	SimulationParameters params;
	dim3 blockSize;
	dim3 numBlocks;
};
