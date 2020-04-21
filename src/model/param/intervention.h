#pragma once

#include <string>
#include <vector>

#include "../cell_model.cuh"

using std::vector;
using std::string;

class Intervention {
public:
	Intervention(string property);
	virtual int getNextIterationNumber() = 0;
	virtual void applyNext(CellModel &model) = 0;
protected:
	void updateParameter(string property, double value, CellModel &model);
	void applyIntervention(string property, double value, CellModel &model);
	string property;
};

class CycleIntervention : public Intervention {
public:
	CycleIntervention(string property, vector<double> values, int interval);
	int getNextIterationNumber() override;
	void applyNext(CellModel &model) override;
private:
	vector<double> values;
	int interval;
	int nextIterationNumber;
	int nextValueIndex;
};

class InterpolateIntervention : public Intervention {
public:
	InterpolateIntervention(string property, vector<double> values, vector<int> times);
	int getNextIterationNumber() override;
	void applyNext(CellModel &model) override;
private:
	vector<double> values;
	vector<int> times;
	int nextIterationNumber;
	int leftValueIndex;
};

// class OscillateIntervention : public Intervention {
// 	OscillateIntervention(string property, string type, vector<double> values, vector<double> timingValues);
// 	int getNextIterationNumber() override;
// 	void applyNext(CellModel &model) override;
// };

class ActivateIntervention : public Intervention {
public:
	ActivateIntervention(string property, double value,
		int iterationNumber);
	int getNextIterationNumber() override;
	void applyNext(CellModel &model) override;
private:
	int iterationNumber;
	double value;
};