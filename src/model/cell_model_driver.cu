#include <iostream>
#include <fstream>
#include <chrono>

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
		SimulationParameters &params) {
	int total_time = TIME_DIFF(t0, CURRENT_TIME);
	auto outputStream = std::ofstream(params.output.runtime.file);
	outputStream << "total_time" << std::endl;
	outputStream << total_time << std::endl;
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

	double totalCellsDuration = 0;
	double totalBigCellsDuration = 0;	
	double totalInteractionsDuration = 0;
	double totalPrepareGrowthDuration = 0;
	double totalGrowthInteractionsDuration = 0;
	double totalEnvironmentDuration = 0;
	
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

		totalCellsDuration += model.durationCells.count();
		totalBigCellsDuration += model.durationBigCells.count();	
		totalInteractionsDuration += model.durationInteractions.count();
		totalPrepareGrowthDuration += model.durationPrepareGrowth.count();
		totalGrowthInteractionsDuration += model.durationGrowthInteractions.count();
		totalEnvironmentDuration += model.durationEnvironment.count();

	}

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

	std::cout << "Total runtime: " << duration.count() << std::endl;

	std::cout << "Runtime of cell function: " << totalCellsDuration << std::endl;
	std::cout << "Runtime of big cell function: " << totalBigCellsDuration << std::endl;
	std::cout << "Runtime of interactions function: " << totalInteractionsDuration << std::endl;
	std::cout << "Runtime of prepare growth function: " << totalPrepareGrowthDuration << std::endl;
	std::cout << "Runtime of growth interactions function: " << totalGrowthInteractionsDuration << std::endl;
	std::cout << "Runtime of environment function: " << totalEnvironmentDuration << std::endl;

	std::cout << "Total function runtimes: " << totalCellsDuration + totalBigCellsDuration + totalInteractionsDuration + totalPrepareGrowthDuration + totalGrowthInteractionsDuration + totalEnvironmentDuration << std::endl;

	videoOutput.close();
	statisticsOutput.close();

	writeRuntime(t0, params);

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