#include <iostream>
#include <stdio.h>

#include "cell_model.cuh"

// TODO: Const things

int main() {
	CellModelParams params(1000, 1000);
	params.initial_density = 0.5;
	params.survival_threshold = 3;
	params.energy_loss_p = 0.9;
	params.gather_light_energy_p = 0.5;

	CellModel model(params);
	// model.printCells();

	// for (int i = 0; i < 10; i++) {
		model.simulate(1000);
		// model.printCells();
		// std::cout << std::endl << std::endl;
	// }
}