#pragma once

#include <curand_kernel.h>

struct Cell {
	enum precision {
		nbits_occupied = 1,
		nbits_energy = 8,
		nbits_waste = 8,
		nbits_co2 = 8,
		nbits_chem = 8
	};

	bool occupied: 	nbits_occupied;
	int energy: 	nbits_energy;
	int waste: 		nbits_waste;
	int co2: 		nbits_co2;
	int chem: 		nbits_chem;

	// Simulation:
	curandState randState;
};

std::ostream& operator<<(std::ostream &stream, const Cell &cell);