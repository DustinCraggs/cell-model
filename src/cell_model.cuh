#pragma once

#include "cell.h"

struct CudaParams {
	dim3 blockSize = 1024;
	dim3 numBlocks = 10;

	// TODO: Reduction buffer params
};

struct CellModelParams {
	int w, h, d;
	int nCells;
	float initialDensity;
	int survivalThreshold;
	float energyLossProb;
	float gatherLightEnergyProb;

	CudaParams cudaParams;

	CellModelParams(int w, int h, int d) {
		this->w = w;
		this->h = h;
		this->d = d;
		this->nCells = w * h * d;
	}
};

class CellModel {
public:
	CellModel(CellModelParams params);
	void printCells();
	void simulate(int nIterations);
	void synchronizeData();

	// Statistics:
	int numberOfLivingCells();
	double totalEnergy();
private:
	void initialise();

	Cell *cells;
	CellModelParams params;
	dim3 blockSize;
	dim3 numBlocks;
};
