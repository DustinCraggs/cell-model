#include <curand_kernel.h>

#include "growth.cuh"

#include "util.cuh"

// TODO: Death of big cells

namespace growth {
	
	__device__
	int* getMostAdjacentDirection(GridElement *grid, GridElement &element, ModelParameters &params);

	__device__
	int* getRandomDirection(GridElement &element);

	__device__
	void getTotalResourcesOfCell(GridElement *grid, int idx, int *totalEnergy,
		int *totalChem, int *totalDToxin, int *nCells, ModelParameters &params);

	__device__
	void releaseResourcesAndTerminateCell(GridElement *grid, int idx, ModelParameters &params);

	__device__
	void releaseCellResources(GridElement &element);

	__device__
	void redistributeResources(GridElement *grid, int idx, int totalEnergy,
		int totalChem, int totalDToxin, int nCells, ModelParameters &params);

	__device__
	void prepare(GridElement *grid, GridElement &element, ModelParameters &params, int iter) {
		if (!element.cell.alive || element.cell.is_subcell) {
			return;
		}
		// TODO: Distribute chemicals and toxins
		int nCells = 1;
		int totalEnergy = element.cell.energy;
		GridElement *currentCell = &element;
		GridElement::Position pos = currentCell->position;
		int parent_idx = element.cell.parent_idx;
		while (currentCell->cell.has_subcell) {
			Cell::NextCellOffset offset = currentCell->cell.nextCellOffset;
			pos = currentCell->position;
			// printf("OFFSET: %d,%d,%d", offset.x, offset.y, offset.z)
			currentCell = &grid[util::get_idx(pos, offset, params)];
			if (!currentCell->cell.alive) {
				return;
			}
			totalEnergy += currentCell->cell.energy;
			nCells++;
			if (currentCell->cell.parent_idx != parent_idx) {
				printf("CHANGED CELL PARENT: %d!\n", currentCell->cell.alive);
			}
			if (nCells == params.maxCellSize) {
				if (currentCell->cell.has_subcell) {
					printf("EXC!\n");
				}
				return;
			}
		}
		pos = currentCell->position;
		currentCell->cell.nextCellOffset = {0};

		GridElement *subcell = currentCell;
		// printf("PREP TOTAL: %d\n", totalEnergy);
		if (totalEnergy > params.growthThreshold) {
			int *direction = getMostAdjacentDirection(grid, *subcell, params);
			// int *direction = getRandomDirection(*subcell);
			// int direction[3] = {1,0,0};

			int targetIdx =
				util::get_idx(pos.x + direction[0], pos.y + direction[1], pos.z + direction[2], params);
			if (targetIdx == -1) {
				return;
			}
			GridElement *targetCell = &grid[targetIdx];
			if (!targetCell->cell.alive) {
				// if (iter < 5) {
				// 	printf("CAN GROW: %d\n", targetIdx);
				// }
				subcell->cell.nextCellOffset = {direction[0], direction[1], direction[2]};
				targetCell->canGrow = true;
				targetCell->cell.is_subcell = false;
				targetCell->cell.has_subcell = false;
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
		// Reset losing cells' growth directions

		if (element.canGrow && !element.cell.alive) {
			element.canGrow = false;
			int growth_intents[6] = {0};
			int n_intents = 0;
			static int directions[6][3] = DIRECTIONS_6_WAY_ARRAY;
			for (int i = 0; i < 6; i++) {
				int checkIdx = util::get_idx(element.position, directions[i], params);
				if (checkIdx != -1) {
					GridElement *check = &grid[checkIdx];
					if (check->cell.alive) {
						if (util::get_idx(check->position, check->cell.nextCellOffset, params)
								== element.position.idx) {
							growth_intents[n_intents++] = i;
							check->cell.nextCellOffset = {0};
							// break;
						}
					}
				}
			}
			int choice = 0;
			if (n_intents == 0) {
				printf("No intents found %d\n", element.position.idx);
				return;
			}
			if (n_intents > 1) {
			// if (n_intents != 1) {
				// printf("COLLISION %d\n", n_intents);
				double tieBreak = curand_uniform(&element.cell.randState);
				tieBreak = tieBreak == 1.0 ? 0.0 : tieBreak;
				choice = tieBreak * (n_intents);
			}
			int *choiceDirection = directions[growth_intents[choice]];
			int choiceIdx = util::get_idx(
				element.position,
				choiceDirection,
				params
			);
			// if (!grid[choiceIdx].cell.alive) {
			// 	return;
			// }

			grid[choiceIdx].cell.has_subcell = true;
			grid[choiceIdx].cell.nextCellOffset = {
				-choiceDirection[0],
				-choiceDirection[1],
				-choiceDirection[2]
			};

			element.cell.is_subcell = true;
			element.cell.has_subcell = false;
			element.cell.nextCellOffset = {};
			// Temporary initialisation:
			element.cell.alive = true;
			element.cell.energy = 0;
			element.cell.chem = 0;
			element.cell.dToxin = 10;
			element.cell.ndToxin = 10;
			element.cell.parent_idx = grid[choiceIdx].cell.parent_idx;

		}
		element.canGrow = false;
	}

	// TODO: Consolidate death check for single- and multi-unit cells
	__device__
	void checkDeathAndDistributeResources(GridElement *grid, GridElement &element, ModelParameters &params) {
		if (!element.cell.alive || element.cell.is_subcell) {
			return;
		}
		int totalEnergy, totalChem, totalDToxin, nCells;
		getTotalResourcesOfCell(grid, element.position.idx, &totalEnergy,
			&totalChem, &totalDToxin, &nCells, params);

		// if (totalEnergy < params.energySurvivalThreshold
		// 		|| totalChem < params.chemSurvivalThreshold) {
		// 		// || totalDToxin >= params.dToxinDeathThreshold) {
		// 	releaseResourcesAndTerminateCell(grid, element.position.idx, params);
		// } else {
			redistributeResources(grid, element.position.idx, totalEnergy,
				totalChem, totalDToxin, nCells, params);
		// }
	}

	__device__
	void getTotalResourcesOfCell(GridElement *grid, int idx, int *totalEnergy,
			int *totalChem, int *totalDToxin, int *nCells, ModelParameters &params) {
		GridElement *subcell = &grid[idx];
		int energy = subcell->cell.energy;
		int chem = subcell->cell.chem;
		int toxin = subcell->cell.dToxin;
		int cellCount = 1;
		if (subcell->cell.parent_idx != idx) {
			printf("NEQ0   %d:%d!\n", idx, subcell->cell.parent_idx);
		}
		while (subcell->cell.has_subcell) {
			GridElement::Position pos = subcell->position;
			Cell::NextCellOffset offset = subcell->cell.nextCellOffset;
			subcell = &grid[util::get_idx(pos, offset, params)];
			if (!subcell->cell.alive) {
				break;
			}
			if (subcell->cell.parent_idx != idx) {
				printf("NEQ%d   %d:%d!\n", cellCount, idx, subcell->cell.parent_idx);
			}

			energy += subcell->cell.energy;
			chem += subcell->cell.chem;
			toxin += subcell->cell.dToxin;
			cellCount++;

			if (cellCount == params.maxCellSize) {
				break;
			}
		}
		*totalEnergy = energy;
		*totalChem = chem;
		*totalDToxin = toxin;
		*nCells = cellCount;
	}

	__device__
	void releaseResourcesAndTerminateCell(GridElement *grid, int idx, ModelParameters &params) {
		GridElement *subcell = &grid[idx];
		subcell->cell.alive = false;
		releaseCellResources(*subcell);

		if (subcell->cell.parent_idx != idx) {
			printf("NEQ0   %d:%d!\n", idx, subcell->cell.parent_idx);
		}

		int nCells = 1;
		while (subcell->cell.has_subcell) {
			// Move to next cell:
			GridElement::Position pos = subcell->position;
			Cell::NextCellOffset offset = subcell->cell.nextCellOffset;
			subcell = &grid[util::get_idx(pos, offset, params)];
			if (!subcell->cell.alive) {
				break;
			}
			if (subcell->cell.parent_idx != idx) {
				printf("NEQ%d   %d:%d!\n", nCells, idx, subcell->cell.parent_idx);
			}
			subcell->cell.alive = false;
			releaseCellResources(*subcell);

			nCells++;
			if (nCells == params.maxCellSize) {
				if (subcell->cell.has_subcell) {
					pos = subcell->position;
					offset = subcell->cell.nextCellOffset;
					subcell = &grid[util::get_idx(pos, offset, params)];
					printf("EXC DEATH %d, %d, a: %d\n", idx, subcell->cell.parent_idx, subcell->cell.alive);
				}
				// break;
			}
		}
	}

		// GridElement *currentCell = &element;
		// GridElement::Position pos = currentCell->position;
		// int parent_idx = element.cell.parent_idx;
		// while (currentCell->cell.has_subcell) {
		// 	pos = currentCell->position;
		// 	Cell::NextCellOffset offset = currentCell->cell.nextCellOffset;
		// 	// printf("OFFSET: %d,%d,%d", offset.x, offset.y, offset.z)
		// 	currentCell = &grid[util::get_idx(pos, offset, params)];
		// 	if (!currentCell->cell.alive) {
		// 		return;
		// 	}
		// 	totalEnergy += currentCell->cell.energy;
		// 	nCells++;
		// 	if (currentCell->cell.parent_idx != parent_idx) {
		// 		printf("CHANGED CELL PARENT: %d!\n", currentCell->cell.alive);
		// 	}
		// 	if (nCells == params.maxCellSize) {
		// 		if (currentCell->cell.has_subcell) {
		// 			printf("EXC!\n");
		// 		}
		// 		return;
		// 	}
		// }

	__device__
	void releaseCellResources(GridElement &element) {
		// Release 90% of resources to env.:
		int maxChem = (1 << Environment::nbits_chem) - 1;
		int newChem = 0.9 * element.cell.chem + element.environment.chem;
		element.environment.chem = newChem <= maxChem ? newChem : maxChem;

		int maxDToxin = (1 << Environment::nbits_d_toxin) - 1;
		int newDToxin = 0.9 * element.cell.dToxin + element.environment.dToxin;
		element.environment.dToxin = newDToxin <= maxDToxin ? newDToxin : maxDToxin;

		int maxNDToxin = (1 << Environment::nbits_nd_toxin) - 1;
		int newNDToxin = 0.9 * element.cell.ndToxin + element.environment.ndToxin;
		element.environment.ndToxin = newNDToxin <= maxNDToxin ? newNDToxin : maxNDToxin;
	}

	__device__
	void redistributeResources(GridElement *grid, int idx, int totalEnergy,
			int totalChem, int totalDToxin, int nCells, ModelParameters &params) {
		GridElement *subcell = &grid[idx];
		int newEnergy = totalEnergy / nCells;
		int newChem = totalChem / nCells;
		int newDToxin = totalDToxin / nCells;

		// Add any leftover energy to first cell:
		subcell->cell.energy = newEnergy + (totalEnergy % nCells);
		subcell->cell.chem = newChem + (totalChem % nCells);
		subcell->cell.dToxin = newDToxin + (totalDToxin % nCells);
		int cellCount = 1;
		while (subcell->cell.has_subcell) {
			// Move to next cell:
			GridElement::Position pos = subcell->position;
			Cell::NextCellOffset offset = subcell->cell.nextCellOffset;
			subcell = &grid[util::get_idx(pos, offset, params)];
			if (!subcell->cell.alive) {
				break;
			}
			// Redistribute resources & toxins:
			subcell->cell.energy = newEnergy;
			subcell->cell.chem = newChem;
			subcell->cell.dToxin = newDToxin;

			cellCount++;
			if (cellCount == params.maxCellSize) {
				break;
			}
		}
	}

	// __device__
	// void getNextSubcell(GridElement)

	__device__
	int* getRandomDirection(GridElement &element) {
		double directionRoll = curand_uniform(&element.cell.randState);
		directionRoll = directionRoll == 1.0 ? 0.0 : directionRoll;
		int idx = (int) (directionRoll * 6);
		return util::DIRECTIONS_6_WAY[idx];
	}


	__device__
	int* getMostAdjacentDirection(GridElement *grid, GridElement &element, ModelParameters &params) {
		int maxAdjacencyCount = 1;
		int* maxAdjacencyDirection = nullptr;

		for (int *dir : util::DIRECTIONS_6_WAY) {
			int idx = util::get_idx(element.position, dir, params);
			if (idx != -1 && !grid[idx].cell.alive) {
				int adjacencyCount = 0;
				for (int *checkDir : util::DIRECTIONS_6_WAY) {
					int checkIdx = util::get_idx(grid[idx].position, checkDir, params);
					if (grid[checkIdx].cell.alive
							&& grid[checkIdx].cell.parent_idx == element.cell.parent_idx) {
						adjacencyCount++;
					}
				}
				// TODO: Pick a random one if even:
				if (adjacencyCount > maxAdjacencyCount) {
					maxAdjacencyCount = adjacencyCount;
					maxAdjacencyDirection = dir;
				}
			}
		}
		if (maxAdjacencyDirection == nullptr) {
			return getRandomDirection(element);
		}
		return maxAdjacencyDirection;
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