#pragma once

#include "../grid_element.h"
#include "../cell_model_params.h"

__device__
void set_move_direction(GridElement *grid, GridElement &element, int x, int y, int z, CellModelParams &params);

__device__
void move_cells(GridElement *grid, GridElement &element, int x, int y, int z, CellModelParams &params);