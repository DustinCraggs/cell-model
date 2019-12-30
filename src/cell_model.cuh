#pragma once

#include "cell.h"

struct CellModelParams {
	int w, h;
	float initial_density;
	int survival_threshold;
	float energy_loss_p;
	float gather_light_energy_p;
	int nThreads = 32;

	CellModelParams(int w, int h) {
		this->w = w;
		this->h = h;
	}
};

class CellModel {
public:
	CellModel(CellModelParams params);
	void printCells();
	void simulate(int nIterations);
private:
	void initialise();

	Cell *cells;
	CellModelParams params;
	dim3 blockSize;
	dim3 numBlocks;
};
