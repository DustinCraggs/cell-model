#include <curand_kernel.h>

#include <stdio.h>

#include "movement.cuh"
#include "../grid_element.h"

#define MOVEMENT_DIRECTIONS {{1, 0}, {-1, 0}, {0, 1}, {0, -1}}

__device__
double get_move_position(GridElement element, CellModelParams &params, int *x, int *y);

__device__
extern int get_idx(int x, int y, int z, CellModelParams &params);

__device__
void set_move_direction(GridElement *grid, GridElement &element, int x, int y, int z, CellModelParams &params) {
	int new_x = x, new_y = y;
	int tie_breaker = get_move_position(element, params, &new_x, &new_y);
	int new_idx = get_idx(new_x, new_y, z, params);

	// Check if new position valid:
	if ((new_idx < 0) || (new_idx >= params.gridSize)) {
		return;
	}
	// Check if new position not occupied:
	if (grid[new_idx].cell.alive) {
		return;
	}
	// Check if another cell is moving into the same square:
	int directions[4][2] = MOVEMENT_DIRECTIONS;
	for (int *dir : directions) {
		int check_x = new_x + dir[0];
		int check_y = new_y + dir[1];

		if (!(check_x == x && check_y == y)) {
			int check_idx = get_idx(check_x, check_y, z, params);
			if (!grid[check_idx].cell.alive) {
				continue;
			}
			if ((check_idx < 0) || (check_idx >= params.gridSize)) {
				continue;
			}
			
			int check_tie_breaker = get_move_position(grid[check_idx], params, &check_x, &check_y);
			bool isSamePosition = (check_x == new_x) && (check_y == new_y);
			if (isSamePosition && (check_tie_breaker >= tie_breaker)) {
				return;
			}
		}
	}

	// This cell can move:
	curandState tempState = element.randState;
	curand_uniform(&tempState);
	if (curand_uniform(&tempState) < params.movementProbability) {
		element.canMove = true;
	}
}

__device__
void move_cells(GridElement *grid, GridElement &element, int x, int y, int z, CellModelParams &params) {
	if (element.canMove) {
		get_move_position(element, params, &x, &y);

		// Update new cell:
		int new_idx = get_idx(x, y, z, params);
		grid[new_idx].cell = element.cell;

		// Update old cell:
		element.cell = {};
		element.canMove = false;
	}
	for (int i = 0; i < 2; i++) {
		curand_uniform(&element.randState);
	}
}

__device__
double get_move_position(GridElement element, CellModelParams &params, int *x, int *y) {
	curandState tempState = element.randState;
	double rand = curand_uniform(&tempState);
	if (rand < 0.25) {
		*y -= 1;
	} else if (rand < 0.5) {
		*y += 1;
		rand -= 0.25;
	} else if (rand < 0.75) {
		*x -= 1;
		rand -= 0.5;
	} else {
		*x += 1;
		rand -= 0.75;
	}
	return rand;
}
