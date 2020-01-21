#pragma once
#pragma pack(push, 1)

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

	// TOT: 33
	// TODO: Remove tests:
	// int t1: 5;
	// int t2: 5;
	// int t3: 5;
	// int t4: 5;
	// int t6: 5;
	// 58
	// int t7: 5;
	// int t8: 5;
	// 68
	// int t9: 5;
	// int t10: 5;
	// int t11: 5;
	// int t12: 5;
	// 88

	// Simulation:
	curandState randState;
};

#pragma pack(pop)

std::ostream& operator<<(std::ostream &stream, const Cell &cell);