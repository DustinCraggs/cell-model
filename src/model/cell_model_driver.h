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
	void writeRuntime(SimulationParameters &params, std::vector<int> runtimeVector);
};