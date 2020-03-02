#include "util.cuh"

#include "../param/simulation_parameters.h"
#include "../grid_element.h"

__host__ __device__
int util::get_idx(int x, int y, int z, ModelParameters &params) {
	bool inBounds = x >= 0 && x < params.w
				 && y >= 0 && y < params.h
				 && z >= 0 && z < params.d;
	return inBounds ? x + y * params.w + z * params.w * params.h : -1;
}

__host__ __device__
int util::get_idx(int x, int y, int z, int *xyz_offset, ModelParameters &params) {
	return util::get_idx(x + xyz_offset[0], y + xyz_offset[1], z + xyz_offset[2], params);
}

__host__ __device__
int util::get_idx(GridElement::Position p, ModelParameters &params) {
	return util::get_idx(p.x, p.y, p.z, params);
}

__host__ __device__
int util::get_idx(GridElement::Position p, int *xyz_offset, ModelParameters &params) {
	return util::get_idx(p.x + xyz_offset[0], p.y + xyz_offset[1], p.z + xyz_offset[2], params);
}


__host__ __device__
int util::get_idx(GridElement::Position p, Cell::NextCellOffset offset, ModelParameters &params) {
	return util::get_idx(p.x + offset.x, p.y + offset.y, p.z + offset.z, params);
}

__host__ __device__
bool operator==(const GridElement::Position& lhs,
		const GridElement::Position& rhs) {
    return lhs.x == rhs.x
	    && lhs.y == rhs.y
	    && lhs.z == rhs.z
	    && lhs.idx == rhs.idx;
}

__host__ __device__
bool operator!=(const GridElement::Position& lhs,
		const GridElement::Position& rhs) {
    return !(lhs == rhs);
}
