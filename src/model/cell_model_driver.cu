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
	}

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