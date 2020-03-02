#include <curand_kernel.h>

#include <stdio.h>

#include "movement.cuh"

#include "util.cuh"
#include "../grid_element.h"
#include "../param/simulation_parameters.h"

#define N_RANDOM_NUMBERS_USED 3

__device__ double get_intended_move(GridElement element,
	ModelParameters &params, GridElement::Position &newPosition);
__device__ void move_position(int direction[3],
	GridElement::Position &position, ModelParameters &params);

namespace movement {

	__device__
	void prepare(GridElement *grid, GridElement &element, ModelParameters &params) {
		if (element.cell.has_subcell || element.cell.is_subcell) {
			return;
		}
		element.canMove = false;
		GridElement::Position newPosition;
		double tieBreaker = get_intended_move(element, params, newPosition);

		if (tieBreaker == -1) {
			// Cell is dead, or new position is out of bounds, or cell is
			// not moving this iteration 
			return;
		}

		if (grid[newPosition.idx].cell.alive) {
			// New position is already occupied on this iteration
			return;
		}

		// Check if any other adjacent cell intends to move into newPosition
		static int directions[6][3] = DIRECTIONS_6_WAY;
		for (int *direction : directions) {
			GridElement::Position otherPosition = newPosition;
			move_position(direction, otherPosition, params);

			if (otherPosition.idx != -1 && otherPosition != element.position) {
				double otherTieBreaker = get_intended_move(
					grid[otherPosition.idx],
					params,
					otherPosition
				);
				if (otherPosition == newPosition && otherTieBreaker >= tieBreaker) {
					// Another cell will move into this cell's intended
					// position, and the other cell has a higher random tie
					// breaker
					return;
				}
			}
		}
		element.canMove = true;
	}

	__device__
	void execute(GridElement *grid, GridElement &element, ModelParameters &params) {
		if (element.canMove) {
			element.canMove = false;
			GridElement::Position newPosition;
			get_intended_move(element, params, newPosition);
			
			// Update new cell:
			grid[newPosition.idx].cell = element.cell;

			// Update old cell:
			element.cell.alive = false;
		}

		// Update this GridElement's random number generator:
		for (int i = 0; i < N_RANDOM_NUMBERS_USED; i++) {
			curand_uniform(&element.randState);
		}
	}
}

__device__
double get_intended_move(GridElement element, ModelParameters &params, 
		GridElement::Position &newPosition) {
	if (!element.cell.alive) {
		return -1;
	}

	// Copy the GridElement's random state to avoid changing it:
	curandState tempState = element.randState;

	double moveRoll = curand_uniform(&tempState);
	if (moveRoll >= params.movementProbability) {
		// No movement:
		return -1;
	}

	double directionRoll = curand_uniform(&tempState);
	if (directionRoll == 1.0) {
		directionRoll = 0.0;
	}
	static int directions[6][3] = DIRECTIONS_6_WAY;
	int *direction = directions[(int) (directionRoll * 7)];
	newPosition = element.position;
	move_position(direction, newPosition, params);
	
	if (newPosition.idx == -1) {
		// Position out of bounds, no movement:
		return -1;
	}

	double tieBreakRoll = curand_uniform(&tempState);
	return tieBreakRoll;
}

__device__
void move_position(int direction[3], GridElement::Position &position,
		ModelParameters &params) {
	position.x = position.x + direction[0];
	position.y = position.y + direction[1];
	position.z = position.z + direction[2];
	position.idx = util::get_idx(position.x, position.y, position.z,
		params);
}
