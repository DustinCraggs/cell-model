#include <iostream>
#include <fstream>
#include <chrono>

#include "cell_model_driver.h"
#include "cell_model.cuh"
#include "param/simulation_parameters.h"
#include "output/video.h"
#include "output/statistics.cuh"

#define CURRENT_TIME std::chrono::steady_clock::now()
#define TIME_DIFF(T0, T1) std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count()

CellModelDriver::CellModelDriver(std::string paramsPath) {
	this->params = SimulationParameters::fromJson(paramsPath);
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
	for (int i = 0; i < params.model.iterations; i++) {
		// std::cout << "Iteration: " << i << std::endl;
		// TODO: Optimise iteration:
		// Video output:
		if ((i % output.video.interval) == 0) {
			videoOutput.write(model, i);
		}

		// Statistics output:
		if ((i % output.statistics.interval) == 0) {
			statisticsOutput.write(model, i);
		}

		model.simulate(1);
	}

	videoOutput.close();
	statisticsOutput.close();

	int total_time = TIME_DIFF(t0, CURRENT_TIME);
	auto outputStream = std::ofstream(params.output.runtime.file);
	outputStream << "total_time" << std::endl;
	outputStream << total_time << std::endl;
	outputStream.close();

	// std::cout << params.output.runtime.file << std::endl;
	// std::cout << "total_time" << std::endl;
	// std::cout << total_time << std::endl;

	std::cout << "COMPLETE" << std::endl;
}