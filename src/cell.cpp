#pragma once

#include <curand_kernel.h>

struct Cell {
	enum precision {
		nbits_occupied = 1,
		nbits_energy = 6,
		nbits_waste = 6
	};

	bool occupied: nbits_occupied;
	int energy: nbits_energy;
	int waste: nbits_waste;
	curandState randState;
};

std::ostream& operator<<(std::ostream &stream, const Cell &cell);