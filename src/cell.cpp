#include <iostream>

#include "Cell.cuh"

std::ostream& operator<<(std::ostream &stream, const Cell &cell) {
	return stream << "O" << cell.occupied << " E" << cell.energy << " W" << cell.waste;
}