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

	int genome1NumberOfLivingCells(CellModel model);
	int genome2NumberOfLivingCells(CellModel model);
	int genome3NumberOfLivingCells(CellModel model);
	int genome4NumberOfLivingCells(CellModel model);
	int genome5NumberOfLivingCells(CellModel model);
	int genome6NumberOfLivingCells(CellModel model);
	int genome7NumberOfLivingCells(CellModel model);
	int genome8NumberOfLivingCells(CellModel model);
	int genome9NumberOfLivingCells(CellModel model);
	int genome10NumberOfLivingCells(CellModel model);

	double genome1TotalCellEnergy(CellModel model);
	double genome2TotalCellEnergy(CellModel model);
	double genome3TotalCellEnergy(CellModel model);
	double genome4TotalCellEnergy(CellModel model);
	double genome5TotalCellEnergy(CellModel model);
	double genome6TotalCellEnergy(CellModel model);
	double genome7TotalCellEnergy(CellModel model);
	double genome8TotalCellEnergy(CellModel model);
	double genome9TotalCellEnergy(CellModel model);
	double genome10TotalCellEnergy(CellModel model);

	double genome1TotalCellChem(CellModel model);
	double genome2TotalCellChem(CellModel model);
	double genome3TotalCellChem(CellModel model);
	double genome4TotalCellChem(CellModel model);
	double genome5TotalCellChem(CellModel model);
	double genome6TotalCellChem(CellModel model);
	double genome7TotalCellChem(CellModel model);
	double genome8TotalCellChem(CellModel model);
	double genome9TotalCellChem(CellModel model);
	double genome10TotalCellChem(CellModel model);

	double genome1TotalCellToxin(CellModel model);
	double genome2TotalCellToxin(CellModel model);
	double genome3TotalCellToxin(CellModel model);
	double genome4TotalCellToxin(CellModel model);
	double genome5TotalCellToxin(CellModel model);
	double genome6TotalCellToxin(CellModel model);
	double genome7TotalCellToxin(CellModel model);
	double genome8TotalCellToxin(CellModel model);
	double genome9TotalCellToxin(CellModel model);
	double genome10TotalCellToxin(CellModel model);
	

private:
	std::ofstream outputStream;
};