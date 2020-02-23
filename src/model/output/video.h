#include <iostream>
#include <stdlib.h>

#include "../cell_model_params.h"
#include "../cell_model.cuh"

class VideoOutput {
public:
	VideoOutput(CellModelParams params, OutputParams outputParams);
	void writeFrame(CellModel model);
	void close();
private:
	FILE* openVideoOutputProcess(std::string path, int w, int h);
	void getEnergyFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx);
	void getChemFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx);
	void getToxinFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx);

	FILE *energyPipe = nullptr;
	FILE *chemPipe = nullptr;
	FILE *toxinPipe = nullptr;
	int nPixels;
	unsigned char* frameBuffer;
};