#include "intervention.h"

#include <iostream>
#include <string>
#include <vector>

#include "simulation_parameters.h"
#include "../cell_model.cuh"

using std::vector;
using std::string;

Intervention::Intervention(string property) : property(property) {}

void Intervention::applyIntervention(string property, double value, CellModel &model) {
	if (property.compare("invertedChemRedistribution") == 0) {
		model.redistributeChemicals(value, true);
	} else if (property.compare("chemRedistribution") == 0) {
		model.redistributeChemicals(value, false);
	} else {
		updateParameter(property, value, model);
	}
}

void Intervention::updateParameter(string property, double value, CellModel &model) {
	SimulationParameters params = model.getParams();
	if (property.compare("lightIntensity") == 0) {
		std::cout << "UPDATING LIGHTINTENSITY: " << value << std::endl;
		params.model.lightIntensity = value;
	} else if (property.compare("temperature") == 0) {
		params.model.temperature = value;
	} else if (property.compare("energyUsageRate") == 0) {
		params.model.energyUsageRate = value;
	}

	model.updateParams(params);
}

// CYCLE:

CycleIntervention::CycleIntervention(string property, vector<double> values,
	int interval) :
		Intervention(property),
		interval(interval),
		values(values),
		nextIterationNumber(0),
		nextValueIndex(0) {}

int CycleIntervention::getNextIterationNumber() {
	return nextIterationNumber;
}

void CycleIntervention::applyNext(CellModel &model) {
	if (nextValueIndex < values.size()) {
		applyIntervention(property, values.at(nextValueIndex++), model);
		nextIterationNumber += interval;
	}

	if (nextValueIndex == values.size()) {
		nextValueIndex = 0;
	}
}

// INTERPOLATE:


InterpolateIntervention::InterpolateIntervention(string property, vector<double> values,
	vector<int> times) :
		Intervention(property),
		values(values),
		times(times),
		nextIterationNumber(0),
		leftValueIndex(0) {
	if (values.size() != times.size()) {
		std::cerr <<
			"ERROR: InterpolateIntervention values and times must be same size"
			<< std::endl;
		std::exit(EXIT_FAILURE);
	}
}

int InterpolateIntervention::getNextIterationNumber() {
	return nextIterationNumber;
}

void InterpolateIntervention::applyNext(CellModel &model) {
	if (leftValueIndex == times.size() - 1) {
		// No more values/times to interpolate
		return;
	}

	// Compute current interpolation:
	double iterationsSinceLeftValue = nextIterationNumber - times.at(leftValueIndex);
	double iterationSpan = times.at(leftValueIndex + 1) - times.at(leftValueIndex);

	double fraction = iterationsSinceLeftValue/iterationSpan;
	double interpolatedValue = values.at(leftValueIndex + 1) * fraction
		+ values.at(leftValueIndex) * (1 - fraction);
	// std::cout << interpolatedValue << std::endl;
	applyIntervention(property, interpolatedValue, model);

	// Update 
	nextIterationNumber++;
	if (times.at(leftValueIndex + 1) < nextIterationNumber) {
		leftValueIndex++;
	}
}

// OSCILLATE:

// OscillateIntervention::OscillateIntervention(string property, string type,
// 		vector<double> values) : Intervention(property, type, values) {

// }

// int OscillateIntervention::getNextIterationNumber() {
// 	return 0;
// }

// void OscillateIntervention::applyNext(SimulationParameters &params, CellModel model) {

// }

// ACTIVATE:

ActivateIntervention::ActivateIntervention(string property, double value,
	int iterationNumber) :
		Intervention(property),
		value(value),
		iterationNumber(iterationNumber) {

}

int ActivateIntervention::getNextIterationNumber() {
	return iterationNumber;
}

void ActivateIntervention::applyNext(CellModel &model) {
	applyIntervention(property, value, model);
	iterationNumber = -1;
}