#pragma once

#include "../grid_element.h"
#include "../param/simulation_parameters.h"

namespace movement {
	__device__
	void prepare(GridElement *grid, GridElement &element, ModelParameters &params);

	__device__
	void execute(GridElement *grid, GridElement &element, ModelParameters &params);
}