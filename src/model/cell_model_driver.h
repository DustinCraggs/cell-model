#pragma once
#include <string>
#include <chrono>
#include <vector>

#include "param/simulation_parameters.h"

class CellModelDriver {
public:
	CellModelDriver(std::string paramsPath);
	void run();
private:
	SimulationParameters params;
	void writeRuntime(std::chrono::time_point<std::chrono::steady_clock> t0,
		SimulationParameters &params, std::vector<double> runtimeVector);
};