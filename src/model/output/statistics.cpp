// #include <thrust/transform_reduce.h>
// #include <thrust/functional.h>

// struct CellOccupied {
// 	__device__
// 	int operator()(const GridElement& g) const {
// 		return g.cell.alive ? 1 : 0;
// 	}
// };

// struct CellEnergy {
// 	__device__
// 	double operator()(const GridElement& g) const {
// 		return g.cell.alive ? g.cell.energy : 0.0;
// 	}
// };

// struct CellChem {
// 	__device__
// 	double operator()(const GridElement& g) const {
// 		return g.cell.alive ? g.cell.chem : 0.0;
// 	}
// };

// struct CellToxin {
// 	__device__
// 	double operator()(const GridElement& g) const {
// 		return g.cell.alive ? g.cell.dToxin + g.cell.ndToxin : 0.0;
// 	}
// };

// struct EnvChem {
// 	__device__
// 	double operator()(const GridElement& g) const {
// 		return g.environment.chem;
// 	}
// };

// struct EnvToxin {
// 	__device__
// 	double operator()(const GridElement& g) const {
// 		return g.environment.dToxin + g.environment.ndToxin;
// 	}
// };

// int CellModel::numberOfLivingCells() {
// 	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
// 	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellOccupied(), 0, thrust::plus<int>());
// }

// double CellModel::totalCellEnergy() {
// 	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
// 	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellEnergy(), 0.0, thrust::plus<double>());
// }

// double CellModel::totalCellChem() {
// 	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
// 	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellChem(), 0.0, thrust::plus<double>());
// }

// double CellModel::totalCellToxin() {
// 	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
// 	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, CellToxin(), 0.0, thrust::plus<double>());
// }

// double CellModel::totalEnvChem() {
// 	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
// 	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, EnvChem(), 0.0, thrust::plus<double>());
// }

// double CellModel::totalEnvToxin() {
// 	thrust::device_ptr<GridElement> gridPtr = thrust::device_pointer_cast(grid);
// 	return thrust::transform_reduce(gridPtr, gridPtr + params.gridSize, EnvToxin(), 0.0, thrust::plus<double>());
// }
