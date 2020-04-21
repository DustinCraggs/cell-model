#pragma once

#include "../param/simulation_parameters.h"
#include "../grid_element.h"

#define DIRECTIONS_6_WAY_ARRAY {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}

namespace util {
	__device__
	static int DIRECTIONS_6_WAY[6][3] = DIRECTIONS_6_WAY_ARRAY;
	
	__host__ __device__
	int get_idx(int x, int y, int z, ModelParameters &params);

	__host__ __device__
	int get_idx(int x, int y, int z, int *xyz_offset, ModelParameters &params);

	__host__ __device__
	int get_idx(GridElement::Position p, ModelParameters &params);

	__host__ __device__
	int get_idx(GridElement::Position p, int *xyz_offset, ModelParameters &params);
	
	__host__ __device__
	int get_idx(GridElement::Position p, Cell::NextCellOffset offset, ModelParameters &params);
}

__host__ __device__
bool operator==(const GridElement::Position& lhs,
		const GridElement::Position& rhs);

__host__ __device__
bool operator!=(const GridElement::Position& lhs,
		const GridElement::Position& rhs);

