#pragma once
#include <string>

#include "param/simulation_parameters.h"

class CellModelDriver {
public:
	CellModelDriver(std::string paramsPath);
	void run();
private:
	SimulationParameters params;
};