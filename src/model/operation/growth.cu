#include <curand_kernel.h>

#include "growth.cuh"

#include "util.cuh"

// TODO: Death of big cells

namespace growth {
	
	__device__
	void prepare(GridElement *grid, GridElement &element, ModelParameters &params) {
		if (element.cell.is_subcell) {
			return;
		}
		// TODO: Distribute chemicals and toxins
		int nCells = 1;
		int totalEnergy = element.cell.energy;
		GridElement &subcell = element;
		GridElement::Position pos = subcell.position;
		while (subcell.cell.has_subcell) {
			Cell::NextCellOffset offset = subcell.cell.nextCellOffset;
			pos = subcell.position;
			// printf("OFFSET: %d,%d,%d", offset.x, offset.y, offset.z)
			subcell = grid[util::get_idx(pos, offset, params)];
			totalEnergy += subcell.cell.energy;
			nCells++;
			if (nCells == params.maxCellSize) {
				return;
			}
		}
		pos = subcell.position;

		// printf("PREP TOTAL: %d\n", totalEnergy);
		if (totalEnergy > params.growthThreshold) {
			static int directions[6][3] = DIRECTIONS_6_WAY;
			double directionRoll = curand_uniform(&subcell.cell.randState);
			directionRoll = directionRoll == 1.0 ? 0.0 : directionRoll;
			int *direction = directions[(int) (directionRoll * 7)];
			int targetIdx =
				util::get_idx(pos.x + direction[0], pos.y + direction[1], pos.z + direction[2], params);
			if (targetIdx == -1) {
				return;
			}
			GridElement &targetCell = 
				grid[targetIdx];
			if (!targetCell.cell.alive) {
				// printf("CAN GROW: %d\n", targetIdx);
				subcell.cell.nextCellOffset = {direction[0], direction[1], direction[2]};
				targetCell.canGrow = true;
			}
		}
	}

	__device__
	void execute(GridElement *grid, GridElement &element, ModelParameters &params) {
		// Execute:
		// Grow empty cell:
		// Check canGrow
		// Check adjacent cells
		// If multiple cells intend to grow into this cell, roll tie break on this empty cell
		// Grow winning cell here (empty until next growth)
		// Reset losing cells' growth direction

		if (element.canGrow) {
			int growth_intents[6] = {0};
			int n_intents = 0;
			static int directions[6][3] = DIRECTIONS_6_WAY;
			for (int i = 0; i < 6; i++) {
				int *offset = directions[i];
				int checkIdx = util::get_idx(element.position, offset, params);
				if (checkIdx != -1) {
					GridElement &check = grid[checkIdx];
					if (util::get_idx(check.position, check.cell.nextCellOffset, params)
							== element.position.idx) {
						growth_intents[n_intents++] = i;
						check.cell.nextCellOffset = {};
					}
				}
			}
			int choice = 0;
			if (n_intents > 1) {
			// if (n_intents != 1) {
				printf("COLLISION %d\n", n_intents);
				double tieBreak = curand_uniform(&element.cell.randState);
				tieBreak = tieBreak == 1.0 ? 0.0 : tieBreak;
				choice = tieBreak * (n_intents + 1);
			}
			int *choiceDirection = directions[growth_intents[choice]];
			int choiceIdx = util::get_idx(
				element.position,
				choiceDirection,
				params
			);
			grid[choiceIdx].cell.has_subcell = true;
			grid[choiceIdx].cell.nextCellOffset = {
				-choiceDirection[0],
				-choiceDirection[1],
				-choiceDirection[2]
			};

			element.cell.is_subcell = true;
			// Temporary initialisation:
			element.cell.alive = true;
			element.cell.energy = 0;
			element.cell.chem = 0;
			element.cell.dToxin = 10;
			element.cell.ndToxin = 10;
			element.cell.parent_idx = grid[choiceIdx].cell.parent_idx;
			// printf("INITIALISING %d\n", choiceIdx);

			element.canGrow = false;

			// TODO: Distribute chemicals and toxins
			int nCells = 1;
			GridElement &subcell = grid[element.cell.parent_idx];
			int totalEnergy = subcell.cell.energy;
			int totalChem = subcell.cell.chem;
			GridElement::Position pos = subcell.position;
			while (subcell.cell.has_subcell) {
				Cell::NextCellOffset offset = subcell.cell.nextCellOffset;
				pos = subcell.position;
				subcell = grid[util::get_idx(pos, offset, params)];

				totalEnergy += subcell.cell.energy;
				totalChem += subcell.cell.chem;
				nCells++;
				if (nCells == params.maxCellSize) {
					break;
				}
			}

			totalEnergy -= params.growthCost;

			int energyPerCell = totalEnergy / (float) nCells;
			int chemPerCell = totalChem / (float) nCells;

			subcell = grid[element.cell.parent_idx];
			pos = subcell.position;
			subcell.cell.energy = energyPerCell;
			subcell.cell.chem = chemPerCell;
			nCells = 1;
			while (subcell.cell.has_subcell) {
				Cell::NextCellOffset offset = subcell.cell.nextCellOffset;
				pos = subcell.position;
				subcell = grid[util::get_idx(pos, offset, params)];

				subcell.cell.energy = energyPerCell;
				subcell.cell.chem = chemPerCell;
				nCells++;
				if (nCells == params.maxCellSize) {
					break;
				}
			}

		}
	}

}

// if (element.cell.energy < params.energySurvivalThreshold
// 		|| element.cell.chem < params.chemSurvivalThreshold
// 		|| element.cell.dToxin >= params.dToxinDeathThreshold) {
// 	element.cell.alive = false;
// 	// Release 90% of resources to env.:
// 	int maxChem = (1 << Environment::nbits_chem) - 1;
// 	int newChem = 0.9 * element.cell.chem + element.environment.chem;
// 	element.environment.chem = newChem <= maxChem ? newChem : maxChem;

// 	int maxDToxin = (1 << Environment::nbits_d_toxin) - 1;
// 	int newDToxin = 0.9 * element.cell.dToxin + element.environment.dToxin;
// 	element.environment.dToxin = newDToxin <= maxDToxin ? newDToxin : maxDToxin;

// 	int maxNDToxin = (1 << Environment::nbits_nd_toxin) - 1;
// 	int newNDToxin = 0.9 * element.cell.ndToxin + element.environment.ndToxin;
// 	element.environment.ndToxin = newNDToxin <= maxNDToxin ? newNDToxin : maxNDToxin;

// 	// TODO: Waste, energy
// }

// Iterate over subcells:
		// Accumulate energy
		// Accumulate shape
		// If energy > threshold
		// 		deduct energy cost from total energy
		// 		pick next growth element (unoccupied, maintain rectangular)
		// 		flag element (will be created in execute)
		// 		redistribute energy evenly