#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_ptr.h>
#include <iostream>
#include <fstream>

#include "statistics.cuh"
#include "../cell_model.cuh"
#include "../grid_element.h"
#include "../param/simulation_parameters.h"
#include <string>

StatisticsOutput::StatisticsOutput(SimulationParameters params) {
	outputStream = std::ofstream(params.output.statistics.file);
	// Write header:

	std::string output = "iteration, number_of_cells, average_cell_size, average_cell_energy, average_cell_chem, average_cell_toxin, total_environment_chem, total_environment_toxin";

	int genomeNum = 10;
	std::string genomeString;

	for(int i = 1; i <= genomeNum; i++) {

		output += ", genome";
		output += std::to_string(i);
		output += "numberOfLivingCells";

	}

	for(int i = 1; i <= genomeNum; i++) {

		output += ", genome";
		output += std::to_string(i);
		output += "AverageCellEnergy";

	}

	for(int i = 1; i <= genomeNum; i++) {

		output += ", genome";
		output += std::to_string(i);
		output += "AverageCellChem";

	}

	outputStream << output << std::endl;


}

void StatisticsOutput::write(CellModel model, int iteration) {

	double n_living_cells = numberOfLivingCells(model);

	std::string output = 	std::to_string(iteration) + ',' + 
							std::to_string(n_living_cells) + ',' + 
							std::to_string(numberOfOccupiedElements(model)/(double)n_living_cells) + ',' +
							std::to_string(totalCellEnergy(model)/n_living_cells) + ',' +
							std::to_string(totalCellChem(model)/n_living_cells) + ',' +
							std::to_string(totalCellToxin(model)/n_living_cells) + ',' +
							std::to_string(totalEnvChem(model)) + ',' +
							std::to_string(totalEnvToxin(model));

	int genome1Cells = genome1NumberOfLivingCells(model);
	int genome2Cells = genome2NumberOfLivingCells(model);
	int genome3Cells = genome3NumberOfLivingCells(model);
	int genome4Cells = genome4NumberOfLivingCells(model);
	int genome5Cells = genome5NumberOfLivingCells(model);
	int genome6Cells = genome6NumberOfLivingCells(model);
	int genome7Cells = genome7NumberOfLivingCells(model);
	int genome8Cells = genome8NumberOfLivingCells(model);
	int genome9Cells = genome9NumberOfLivingCells(model);
	int genome10Cells = genome10NumberOfLivingCells(model);

	output += ',' + std::to_string(genome1Cells);
	output += ',' + std::to_string(genome2Cells);
	output += ',' + std::to_string(genome3Cells);
	output += ',' + std::to_string(genome4Cells);
	output += ',' + std::to_string(genome5Cells);
	output += ',' + std::to_string(genome6Cells);
	output += ',' + std::to_string(genome7Cells);
	output += ',' + std::to_string(genome8Cells);
	output += ',' + std::to_string(genome9Cells);
	output += ',' + std::to_string(genome10Cells);

	output += ',' + std::to_string(genome1TotalCellEnergy(model)/genome1Cells);
	output += ',' + std::to_string(genome2TotalCellEnergy(model)/genome2Cells);
	output += ',' + std::to_string(genome3TotalCellEnergy(model)/genome3Cells);
	output += ',' + std::to_string(genome4TotalCellEnergy(model)/genome4Cells);
	output += ',' + std::to_string(genome5TotalCellEnergy(model)/genome5Cells);
	output += ',' + std::to_string(genome6TotalCellEnergy(model)/genome6Cells);
	output += ',' + std::to_string(genome7TotalCellEnergy(model)/genome7Cells);
	output += ',' + std::to_string(genome8TotalCellEnergy(model)/genome8Cells);
	output += ',' + std::to_string(genome9TotalCellEnergy(model)/genome9Cells);
	output += ',' + std::to_string(genome10TotalCellEnergy(model)/genome10Cells);

	output += ',' + std::to_string(genome1TotalCellChem(model)/genome1Cells);
	output += ',' + std::to_string(genome2TotalCellChem(model)/genome2Cells);
	output += ',' + std::to_string(genome3TotalCellChem(model)/genome3Cells);
	output += ',' + std::to_string(genome4TotalCellChem(model)/genome4Cells);
	output += ',' + std::to_string(genome5TotalCellChem(model)/genome5Cells);
	output += ',' + std::to_string(genome6TotalCellChem(model)/genome6Cells);
	output += ',' + std::to_string(genome7TotalCellChem(model)/genome7Cells);
	output += ',' + std::to_string(genome8TotalCellChem(model)/genome8Cells);
	output += ',' + std::to_string(genome9TotalCellChem(model)/genome9Cells);
	output += ',' + std::to_string(genome10TotalCellChem(model)/genome10Cells);
	
	// }

	// else {

	// 	return;

	// }

	outputStream << output << std::endl;

}

void StatisticsOutput::close() {
	outputStream.close();
}

struct cellGenome1 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 1) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome2 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 2) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome3 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 3) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome4 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 4) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome5 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 5) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome6 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 6) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome7 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 7) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome8 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 8) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome9 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 9) {
			return 1;
		}
		return 0;
	}
};

struct cellGenome10 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && !g.cell.is_subcell && g.cell.genome == 10) {
			return 1;
		}
		return 0;
	}
};

struct CellEnergyGenome1 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 1) {
			return g.cell.energy;
		}
		return 0.0;
	}
};

struct CellEnergyGenome2 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 2) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome3 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 3) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome4 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 4) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome5 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 5) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome6 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 6) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome7 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 7) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome8 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 8) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome9 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 9) {
			return g.cell.energy;
		}
		return 0.0;
	}
};


struct CellEnergyGenome10 {
	__device__
	double operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 10) {
			return g.cell.energy;
		}
		return 0.0;
	}
};

struct CellChemGenome1 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 1) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome2 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 2) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome3 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 3) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome4 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 4) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome5 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 5) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome6 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 6) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome7 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 7) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome8 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 8) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome9 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 9) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

struct CellChemGenome10 {
	__device__
	int operator()(const GridElement& g) const {
		if(g.cell.alive && g.cell.genome == 10) {
			return g.cell.chem;
		}
		return 0.0;
	}
};

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

int StatisticsOutput::genome1NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome1(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome2NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome2(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome3NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome3(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome4NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome4(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome5NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome5(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome6NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome6(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome7NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome7(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome8NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome8(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome9NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome9(),
		0,
		thrust::plus<int>()
	);
}

int StatisticsOutput::genome10NumberOfLivingCells(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());

	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		cellGenome10(),
		0,
		thrust::plus<int>()
	);
}

double StatisticsOutput::genome1TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome1(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome2TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome2(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome3TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome3(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome4TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome4(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome5TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome5(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome6TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome6(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome7TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome7(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome8TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome8(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome9TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome9(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome10TotalCellEnergy(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellEnergyGenome10(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome1TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome1(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome2TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome2(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome3TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome3(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome4TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome4(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome5TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome5(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome6TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome6(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome7TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome7(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome8TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome8(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome9TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome9(),
		0.0,
		thrust::plus<double>()
	);
}

double StatisticsOutput::genome10TotalCellChem(CellModel model) {
	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(model.getDeviceGrid());
	return thrust::transform_reduce(
		gridPtr,
		gridPtr + model.getParams().gridSize(),
		CellChemGenome10(),
		0.0,
		thrust::plus<double>()
	);
}




