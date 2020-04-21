#pragma once

#include "../grid_element.h"
#include "../param/simulation_parameters.h"


namespace growth {

	__device__
	void prepare(GridElement *grid, GridElement &element, ModelParameters &params, int iter);

	__device__
	void execute(GridElement *grid, GridElement &element, ModelParameters &params);

	__device__
	void checkDeathAndDistributeResources(GridElement *grid, GridElement &element, ModelParameters &params);

}