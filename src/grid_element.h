#pragma once

#include <curand_kernel.h>

struct Cell {
	enum precision {
		// nbits_link_map = 4,
		nbits_alive = 1,
		nbits_energy = 8,
		nbits_waste = 8,
		nbits_co2 = 8,
		nbits_chem = 8,
		nbits_d_toxin = 8,
		nbits_nd_toxin = 8,
		nbits_d_toxin_storage = 8,
		nbits_nd_toxin_storage = 8
	};

	// unsigned int linkMap:	nbits_link_map;

	bool alive: 			nbits_alive;
	unsigned int energy: 	nbits_energy;
	unsigned int co2: 		nbits_co2;
	unsigned int chem: 		nbits_chem;
	unsigned int dToxin: 	nbits_d_toxin;
	unsigned int ndToxin: 	nbits_nd_toxin;

	unsigned int dToxin_storage: 	nbits_d_toxin_storage;
	unsigned int ndToxin_storage: 	nbits_nd_toxin_storage;
};

struct Environment {
	enum precision {
		nbits_chem = 8,
		nbits_d_toxin = 8,
		nbits_nd_toxin = 8
	};

	unsigned int chem: nbits_chem;
	unsigned int dToxin: nbits_d_toxin;
	unsigned int ndToxin: nbits_nd_toxin;
};

struct GridElement {
	Cell cell;
	Environment environment;

	// Simulation state:
	curandState randState;
	bool canMove, canConsumeChem;
};
