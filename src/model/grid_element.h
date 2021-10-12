#pragma once

#include <curand_kernel.h>

struct __align__(16) Cell {
	enum precision {
		nbits_alive = 8,
		nbits_energy = 8,
		nbits_chem = 8,
		nbits_d_toxin = 8,
		nbits_nd_toxin = 8,
		nbits_d_toxin_storage = 8,
		nbits_nd_toxin_storage = 8
	};

	bool alive;
	unsigned int energy;
	unsigned int chem;
	unsigned int dToxin;
	unsigned int ndToxin;
	unsigned int genome;
	unsigned int energyUsageRate;
	unsigned int chemUsageRate;

	unsigned int dToxin_storage;
	unsigned int ndToxin_storage;

	// Multicell:
	bool is_subcell = false;
	int has_subcell = false;
	// TODO: Replace with idx?
	struct NextCellOffset {
		int x, y, z;
	} nextCellOffset;
	int parent_idx = -1;

	// Cell's internal rng (to be used only when operating on this cell)
	curandState randState;
};

struct Environment {
	enum precision {
		nbits_chem = 8,
		nbits_d_toxin = 8,
		nbits_nd_toxin = 8
	};

	unsigned int chem;
	unsigned int dToxin;
	unsigned int ndToxin;

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
	bool canMove = false, canConsumeChem = false, canGrow = false;
};
