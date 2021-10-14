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

#define CURRENT_TIME std::chrono::steady_clock::now()
#define TIME_DIFF(T0, T1) std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count()

void performInterventions(SimulationParameters &params, CellModel &model, int iteration);

CellModelDriver::CellModelDriver(std::string paramsPath) {
	this->params = SimulationParameters::fromJson(paramsPath);
}

void CellModelDriver::writeRuntime(std::chrono::time_point<std::chrono::steady_clock> t0,
		SimulationParameters &params, std::vector<double> runtimeVector) {
	int total_time = TIME_DIFF(t0, CURRENT_TIME);
	auto outputStream = std::ofstream(params.output.runtime.file);
	outputStream << "total_time, cell_runtime, big_cell_runtime, iteractions_runtime, prepare_growth_runtime, growth_interactions_runtime, environment_runtime, total_function_runtime" << std::endl;
	// outputStream << total_time << "," << model.totalCellsDuration << "," << model.totalBigCellsDuration << "," << model.totalInteractionsDuration << "," << model.totalPrepareGrowthDuration << "," << model.totalGrowthInteractionsDuration << "," << model.totalEnvironmentDuration << std::endl;
	std::string outputValues;

	// for(int i = 0; i < runtimeVector.size(); i++) {

	// 	outputValues += 

	// }

	outputStream.close();
}

void CellModelDriver::run() {

	OutputParameters output = params.output;

	auto t0 = CURRENT_TIME;
	// Create model:
	CellModel model(params);

	// Create output modules:
	VideoOutput videoOutput(params);
	StatisticsOutput statisticsOutput(params);
	
	auto start = std::chrono::high_resolution_clock::now();

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

		model.totalCellsDuration += model.durationCells.count();
		model.totalBigCellsDuration += model.durationBigCells.count();	
		model.totalInteractionsDuration += model.durationInteractions.count();
		model.totalPrepareGrowthDuration += model.durationPrepareGrowth.count();
		model.totalGrowthInteractionsDuration += model.durationGrowthInteractions.count();
		model.totalEnvironmentDuration += model.durationEnvironment.count();

	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Total runtime: " << duration.count() << std::endl;

	std::cout << "Runtime of cell function: " << model.totalCellsDuration << std::endl;
	std::cout << "Runtime of big cell function: " << model.totalBigCellsDuration << std::endl;
	std::cout << "Runtime of interactions function: " << model.totalInteractionsDuration << std::endl;
	std::cout << "Runtime of prepare growth function: " << model.totalPrepareGrowthDuration << std::endl;
	std::cout << "Runtime of growth interactions function: " << model.totalGrowthInteractionsDuration << std::endl;
	std::cout << "Runtime of environment function: " << model.totalEnvironmentDuration << std::endl;

	std::cout << "Total function runtimes: " << model.totalCellsDuration + model.totalBigCellsDuration + model.totalInteractionsDuration + model.totalPrepareGrowthDuration + model.totalGrowthInteractionsDuration + model.totalEnvironmentDuration << std::endl;

	videoOutput.close();
	statisticsOutput.close();

	std::vector<double> runtimeVector;

	runtimeVector.push_back(model.totalCellsDuration);
	runtimeVector.push_back(model.totalBigCellsDuration);
	runtimeVector.push_back(model.totalInteractionsDuration);
	runtimeVector.push_back(model.totalPrepareGrowthDuration);
	runtimeVector.push_back(model.totalGrowthInteractionsDuration);
	runtimeVector.push_back(model.totalEnvironmentDuration);

	writeRuntime(t0, params, runtimeVector);

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