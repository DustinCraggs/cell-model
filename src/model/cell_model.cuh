#pragma once

#include "grid_element.h"
#include "cell_model_params.h"

class CellModel {
public:
	CellModel(CellModelParams params);
	void printCells();
	void simulate(int nIterations);
	void synchronizeData();

	// get host frame buffer
	GridElement* getHostGrid();
	// get device frame buffer
	// ... getDeviceGrid();
private:
	void initialise();
	void iterateRandomMovement();

	GridElement *grid;
	CellModelParams params;
	dim3 blockSize;
	dim3 numBlocks;
};
