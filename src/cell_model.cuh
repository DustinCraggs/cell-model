#pragma once

#include "cell.h"

struct CudaParams {
	int blockSize = 1024;
	int numBlocks = 10;
	
	// TODO: Reduction buffer params
}

struct CellModelParams {
	int w, h, d;
	float initialDensity;
	int survivalThreshold;
	float energyLossProb;
	float gatherLightEnergyP;

	CudaParams cudaParams;

	CellModelParams(int w, int h, int d) {
		this->w = w;
		this->h = h;
		this->d = d;
	}
};

class CellModel {
public:
	CellModel(CellModelParams params);
	void printCells();
	void simulate(int nIterations);
	double* getStatistics();
private:
	void initialise();

	Cell *cells;
	CellModelParams params;
	dim3 blockSize;
	dim3 numBlocks;
};
