#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_ptr.h>
#include <iostream>
#include <fstream>

#include "statistics.cuh"
#include "../cell_model.cuh"
#include "../grid_element.h"
#include "../param/simulation_parameters.h"

StatisticsOutput::StatisticsOutput(SimulationParameters params) {
	outputStream = std::ofstream(params.output.statistics.file);
	// Write header:
	outputStream << "iteration,"
		<< "number_of_cells,average_cell_size,average_cell_energy,average_cell_chem,average_cell_toxin,"
		<< "total_environment_chem,total_environment_toxin" << std::endl;
}

void StatisticsOutput::write(CellModel model, int iteration) {
	double n_living_cells = numberOfLivingCells(model);
	outputStream << iteration << ','
		<< n_living_cells << ','
		<< numberOfOccupiedElements(model)/(double)n_living_cells << ','
		<< totalCellEnergy(model)/n_living_cells << ','
		<< totalCellChem(model)/n_living_cells << ','
		<< totalCellToxin(model)/n_living_cells << ','
		<< totalEnvChem(model) << ','
		<< totalEnvToxin(model) << std::endl;
}

void StatisticsOutput::close() {
	outputStream.close();
}

struct IsMaincell {
	__device__
	int operator()(const GridElement& g) const {
		return g.cell.alive && !g.cell.is_subcell ? 1 : 0;
	}
};

struct CellOccupied {
	__device__
	int operator()(const GridElement& g) const {
		return g.cell.alive ? 1 : 0;
	}
};

struct CellEnergy {
	__device__
	double operator()(const GridElement& g) const {
		return g.cell.alive ? g.cell.energy : 0.0;
	}
};

struct CellChem {
	__device__
	double operator()(const GridElement& g) const {
		return g.cell.alive ? g.cell.chem : 0.0;
	}
};

struct CellToxin {
	__device__
	double operator()(const GridElement& g) const {
		return g.cell.alive ? g.cell.dToxin + g.cell.ndToxin : 0.0;
	}
};

struct EnvChem {
	__device__
	double operator()(const GridElement& g) const {
		return g.environment.chem;
	}
};

struct EnvToxin {
	__device__
	double operator()(const GridElement& g) const {
		return g.environment.dToxin + g.environment.ndToxin;
	}
};

int StatisticsOutput::numberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		IsMaincell(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::numberOfOccupiedElements(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellOccupied(),
		0,
		thrust::plus<int>()
	);
}

double StatisticsOutput::totalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergy(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::totalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChem(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::totalCellToxin(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellToxin(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::totalEnvChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		EnvChem(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::totalEnvToxin(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		EnvToxin(),
		0.0,
		thrust::plus<double>()
	);
}
