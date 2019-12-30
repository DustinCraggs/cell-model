#include <iostream>

#include "cell.h"

std::ostream& operator<<(std::ostream &stream, const Cell &cell) {
	return stream << "O" << cell.occupied << " E" << cell.energy << " W" << cell.waste;
}