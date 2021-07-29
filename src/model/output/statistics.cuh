#pragma once

#include <fstream>

#include "../cell_model.cuh"
#include "../param/simulation_parameters.h"

class StatisticsOutput {
public:
	StatisticsOutput(SimulationParameters params);

	void write(CellModel model, int iteration);
	void close();

	int numberOfLivingCells(CellModel model);
	int numberOfOccupiedElements(CellModel model);
	double totalCellEnergy(CellModel model);
	double totalCellChem(CellModel model);
	double totalCellToxin(CellModel model);
	double totalEnvChem(CellModel model);
	double totalEnvToxin(CellModel model);
private:
	std::ofstream outputStream;
};