#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include "cell_model_driver.h"
#include "cell_model.cuh"
#include "param/simulation_parameters.h"
#include "output/video.h"
#include "output/statistics.cuh"

#include "param/intervention.h"

void performInterventions(SimulationParameters &params, CellModel &model, int iteration);

CellModelDriver::CellModelDriver(std::string paramsPath) {
	this->params = SimulationParameters::fromJson(paramsPath);
}

void CellModelDriver::writeRuntime(SimulationParameters &params, std::vector<int> runtimeVector) {
	auto outputStream = std::ofstream(params.output.runtime.file);

	outputStream << "cell_runtime, big_cell_runtime, iteractions_runtime, prepare_growth_runtime, growth_interactions_runtime, environment_runtime, total_function_runtime" << std::endl;
	std::string outputValues;

	int totalRuntime = 0;

	for(int i = 0; i < runtimeVector.size(); i++) {

		outputValues += "," + std::to_string(runtimeVector.at(i));
		totalRuntime += runtimeVector.at(i);

	}

	outputValues += "," + std::to_string(totalRuntime);

	outputStream << outputValues << std::endl;

	outputStream.close();

}

void CellModelDriver::run() {

	OutputParameters output = params.output;

	int setupRuntime;

	auto start = std::chrono::high_resolution_clock::now();

	// Create model:
	CellModel model(params);

	// Create output modules:
	VideoOutput videoOutput(params);
	StatisticsOutput statisticsOutput(params);
	
	auto stop = std::chrono::high_resolution_clock::now();
	setupRuntime = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

	start = std::chrono::high_resolution_clock::now();

	// Iterate:
	for (int i = 0; i <= params.model.iterations; i++) {

		std::cout << "Iteration: " << i << std::endl;

		// Video output:
		if ((i % output.video.interval) == 0) {
			videoOutput.write(model, i);
		}

		// Statistics output:
		if ((i % output.statistics.interval) == 0) {
			statisticsOutput.write(model, i);
		}

		performInterventions(params, model, i);
		model.simulate(1);

		model.totalCellsRuntime += model.runtimeCells;
		model.totalBigCellsRuntime += model.runtimeBigCells;	
		model.totalInteractionsRuntime += model.runtimeInteractions;
		model.totalPrepareGrowthRuntime += model.runtimePrepareGrowth;
		model.totalGrowthInteractionsRuntime += model.runtimeGrowthInteractions;
		model.totalEnvironmentRuntime += model.runtimeEnvironment;

	}

	stop = std::chrono::high_resolution_clock::now();
	int simulationRuntime = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

	start = std::chrono::high_resolution_clock::now();

	videoOutput.close();
	statisticsOutput.close();

	stop = std::chrono::high_resolution_clock::now();

	int finalisationRuntime = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

	std::cout << std::endl;

	std::cout << "Initialisation runtime: " << setupRuntime << std::endl;
	std::cout << "Simulation runtime: " << simulationRuntime << std::endl;
	std::cout << "Finalisation runtime: " << finalisationRuntime << std::endl;
	std::cout << "Total runtime: " << setupRuntime + simulationRuntime + finalisationRuntime << std::endl;

	std::cout << std::endl;

	std::cout << "Runtime of cell function: " << model.totalCellsRuntime << std::endl;
	std::cout << "Runtime of big cell function: " << model.totalBigCellsRuntime << std::endl;
	std::cout << "Runtime of interactions function: " << model.totalInteractionsRuntime << std::endl;
	std::cout << "Runtime of prepare growth function: " << model.totalPrepareGrowthRuntime << std::endl;
	std::cout << "Runtime of growth interactions function: " << model.totalGrowthInteractionsRuntime << std::endl;
	std::cout << "Runtime of environment function: " << model.totalEnvironmentRuntime << std::endl;

	std::cout << std::endl;

	std::cout << "Runtime total of functions: " << model.totalCellsRuntime + model.totalBigCellsRuntime + model.totalInteractionsRuntime + model.totalPrepareGrowthRuntime + model.totalGrowthInteractionsRuntime + model.totalEnvironmentRuntime << std::endl;

	std::vector<int> runtimeVector;

	runtimeVector.push_back(model.totalCellsRuntime);
	runtimeVector.push_back(model.totalBigCellsRuntime);
	runtimeVector.push_back(model.totalInteractionsRuntime);
	runtimeVector.push_back(model.totalPrepareGrowthRuntime);
	runtimeVector.push_back(model.totalGrowthInteractionsRuntime);
	runtimeVector.push_back(model.totalEnvironmentRuntime);

	writeRuntime(params, runtimeVector);

	std::cout << "COMPLETE" << std::endl;
}

void performInterventions(SimulationParameters &params, CellModel &model, int iteration) {
	for (int i = 0; i < params.nInterventions; i++) {
		Intervention *intervention = params.interventions[i];
		if (intervention->getNextIterationNumber() == iteration) {
			intervention->applyNext(model);
		}
	}
}