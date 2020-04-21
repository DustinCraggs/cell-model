#include <iostream>
#include <stdlib.h>

#include "../param/simulation_parameters.h"
#include "../cell_model.cuh"

class VideoOutput {
public:
	VideoOutput(SimulationParameters params);
	void write(CellModel model, int iteration);
	void close();
private:
	FILE* openVideoOutputProcess(std::string path, int w, int h);
	void writeFrame(CellModel model);
	void getEnergyFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx);
	void getChemFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx);
	void getToxinFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx);

	FILE *energyPipe = nullptr;
	FILE *chemPipe = nullptr;
	FILE *toxinPipe = nullptr;
	int nPixels;
	unsigned char* frameBuffer;
	bool active;
};