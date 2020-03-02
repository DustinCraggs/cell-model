#pragma once

#include <curand_kernel.h>

struct Cell {
	enum precision {
		// nbits_link_map = 4,
		nbits_alive = 1,
		nbits_energy = 8,
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

	// Multicell:
	bool is_subcell;
	int has_subcell;
	// TODO: Replace with idx?
	struct NextCellOffset {
		int x, y, z;
	} nextCellOffset;
	int parent_idx;

	// Cell's internal rng (to be used only when operating on this cell)
	curandState randState;
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

	// Environment's internal rng (to be used only when operating on this
	// environment element)
	curandState randState;
};

struct GridElement {
	Cell cell;
	Environment environment;

	// Position:
	struct Position {
		int x, y, z;
		int idx;
	} position;

	// Simulation state:
	curandState randState;

	// Status: (TODO: move and optimise)
	bool canMove, canConsumeChem, canGrow;
};
