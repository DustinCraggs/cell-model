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

	int modelRuntime = 0;
	int videoOutputRuntime = 0;
	int statisticsOutputRuntime = 0;

	// Iterate:
	for (int i = 0; i <= params.model.iterations; i++) {

		std::cout << "Iteration: " << i << std::endl;

		// Video output:
		if (((i % output.video.interval) == 0) && (output.video.enabled == true)) {
			auto start = std::chrono::high_resolution_clock::now();
			videoOutput.write(model, i);
			auto stop = std::chrono::high_resolution_clock::now();
			videoOutputRuntime += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		}

		// Statistics output:
		if (((i % output.statistics.interval) == 0) && (output.statistics.enabled == true)) {
			auto start = std::chrono::high_resolution_clock::now();
			statisticsOutput.write(model, i);
			auto stop = std::chrono::high_resolution_clock::now();
			statisticsOutputRuntime += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();
		}

		performInterventions(params, model, i);

		auto start = std::chrono::high_resolution_clock::now();
		model.simulate(1);
		auto stop = std::chrono::high_resolution_clock::now();
		modelRuntime += std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

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
	std::cout << "Model runtime: " << modelRuntime << std::endl;
	std::cout << "Video output runtime: " << videoOutputRuntime << std::endl;
	std::cout << "Statistics output runtime: " << statisticsOutputRuntime << std::endl;
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

	std::cout << "Runtime of all functions: " << modelRuntime << std::endl;

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