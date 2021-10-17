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

	outputStream << "setup_runtime, simulation_runtime, video_runtime, statistics_runtime, model_runtime, finalisation_runtime, total_runtime, cell_runtime, big_cell_runtime, iteractions_runtime, prepare_growth_runtime, growth_interactions_runtime, environment_runtime, model_runtime" << std::endl;
	std::string outputValues = std::to_string(runtimeVector.at(0));

	for(int i = 1; i < runtimeVector.size(); i++) {

		outputValues += "," + std::to_string(runtimeVector.at(i));

	}

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

	int totalCellsRuntime = 0;
	int totalBigCellsRuntime = 0;	
	int totalInteractionsRuntime = 0;
	int totalPrepareGrowthRuntime = 0;
	int totalGrowthInteractionsRuntime = 0;
	int totalEnvironmentRuntime = 0;

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

		totalCellsRuntime += model.runtimeCells;
		totalBigCellsRuntime += model.runtimeBigCells;	
		totalInteractionsRuntime += model.runtimeInteractions;
		totalPrepareGrowthRuntime += model.runtimePrepareGrowth;
		totalGrowthInteractionsRuntime += model.runtimeGrowthInteractions;
		totalEnvironmentRuntime += model.runtimeEnvironment;

	}

	stop = std::chrono::high_resolution_clock::now();
	int simulationRuntime = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

	start = std::chrono::high_resolution_clock::now();

	videoOutput.close();
	statisticsOutput.close();

	stop = std::chrono::high_resolution_clock::now();

	int finalisationRuntime = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

	std::cout << std::endl;

	std::cout << "Setup runtime: " << setupRuntime << std::endl;
	std::cout << "Simulation runtime: " << simulationRuntime << std::endl;
	std::cout << "Video output runtime: " << videoOutputRuntime << std::endl;
	std::cout << "Statistics output runtime: " << statisticsOutputRuntime << std::endl;
	std::cout << "Finalisation runtime: " << finalisationRuntime << std::endl;
	int totalRuntime = setupRuntime + simulationRuntime + finalisationRuntime;
	std::cout << "Total runtime: " << totalRuntime << std::endl;

	std::cout << std::endl;

	std::cout << "Runtime of cell function: " << totalCellsRuntime << std::endl;
	std::cout << "Runtime of big cell function: " << totalBigCellsRuntime << std::endl;
	std::cout << "Runtime of interactions function: " << totalInteractionsRuntime << std::endl;
	std::cout << "Runtime of prepare growth function: " << totalPrepareGrowthRuntime << std::endl;
	std::cout << "Runtime of growth interactions function: " << totalGrowthInteractionsRuntime << std::endl;
	std::cout << "Runtime of environment function: " << totalEnvironmentRuntime << std::endl;

	std::cout << std::endl;

	std::cout << "Runtime of all functions: " << modelRuntime << std::endl;

	std::vector<int> runtimeVector;

	runtimeVector.push_back(setupRuntime);
	runtimeVector.push_back(simulationRuntime);
	runtimeVector.push_back(videoOutputRuntime);
	runtimeVector.push_back(statisticsOutputRuntime);
	runtimeVector.push_back(modelRuntime);
	runtimeVector.push_back(finalisationRuntime);
	runtimeVector.push_back(totalRuntime);

	runtimeVector.push_back(totalCellsRuntime);
	runtimeVector.push_back(totalBigCellsRuntime);
	runtimeVector.push_back(totalInteractionsRuntime);
	runtimeVector.push_back(totalPrepareGrowthRuntime);
	runtimeVector.push_back(totalGrowthInteractionsRuntime);
	runtimeVector.push_back(totalEnvironmentRuntime);

	runtimeVector.push_back(modelRuntime);

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