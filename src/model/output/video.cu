#include "video.h"

#include <string>
#include <stdlib.h>

#include "../param/simulation_parameters.h"

using std::string;

VideoOutput::VideoOutput(SimulationParameters params) {
	this->active = params.output.video.active;

	energyPipe = openVideoOutputProcess(
		params.output.video.energyPath,
		params.model.w,
		params.model.h
	);
	chemPipe = openVideoOutputProcess(
		params.output.video.chemPath,
		params.model.w,
		params.model.h
	);
	toxinPipe = openVideoOutputProcess(
		params.output.video.toxinPath,
		params.model.w,
		params.model.h
	);
	genomePipe = openVideoOutputProcess(
		params.output.video.genomePath,
		params.model.w,
		params.model.h
	);

	nPixels = params.model.w * params.model.h;
	frameBuffer = new unsigned char[nPixels * 3];
}

void VideoOutput::write(CellModel model, int iteration) {
	if (!active) {
		return;
	}
	writeFrame(model);
}

void VideoOutput::writeFrame(CellModel model) {
	GridElement *grid;
	if (energyPipe || chemPipe || toxinPipe || genomePipe) { 
		grid = model.getHostGrid();
		if (energyPipe) {
			getEnergyFrame(grid, nPixels, frameBuffer, 0);
			fwrite(frameBuffer, sizeof(unsigned char), nPixels * 3, energyPipe);
		}

		if (chemPipe) {
			getChemFrame(grid, nPixels, frameBuffer, 0);
			fwrite(frameBuffer, sizeof(unsigned char), nPixels * 3, chemPipe);
		}

		if (toxinPipe) {
			getToxinFrame(grid, nPixels, frameBuffer, 0);
			fwrite(frameBuffer, sizeof(unsigned char), nPixels * 3, toxinPipe);
		}

		if(genomePipe) {
			getGenomeFrame(grid, nPixels, frameBuffer, 0);
			fwrite(frameBuffer, sizeof(unsigned char), nPixels * 3, genomePipe);
		}
	}
}

void VideoOutput::close() {
	if (!active) {
		return;
	}
	if (energyPipe && pclose(energyPipe) != 0) {
		std::cerr << "Error: Failed to close energy FFmpeg output process" << std::endl;
	}
	if (chemPipe && pclose(chemPipe) != 0) {
		std::cerr << "Error: Failed to close chem FFmpeg output process" << std::endl;
	}
	if (toxinPipe && pclose(toxinPipe) != 0) {
		std::cerr << "Error: Failed to close toxin FFmpeg output process" << std::endl;
	}
	if (genomePipe && pclose(genomePipe) != 0) {
		std::cerr << "Error: Failed to close genome FFmpeg output process" << std::endl;
	}
}

FILE* VideoOutput::openVideoOutputProcess(std::string path, int w, int h) {
	if (path.empty()) {
		return nullptr;
	}

	const string ffmpegCommand =
		string("ffmpeg -hide_banner -loglevel panic -y -f rawvideo -pixel_format rgb24 -video_size ")
		+ std::to_string(w) + string("x") + std::to_string(h)
		+ string(" -i - -c:v h264 -pix_fmt yuv420p -s 512x512 -sws_flags neighbor ")
		+ path;

	FILE *outputFile = popen(ffmpegCommand.c_str(), "w");
	if (!outputFile) {
		std::cerr << "Error: Could not open FFmpeg process for video output" << std::endl;
		return nullptr;
	}
	return outputFile;
}

void VideoOutput::getEnergyFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx) {
	for (int i = 0; i < nPixels; i++) {
		GridElement element = grid[i + dIdx*nPixels];
		if (element.cell.alive) {
			frameBuffer[i*3] = 0.5 * (255 - element.cell.energy);
			frameBuffer[i*3 + 1] = 0;
			frameBuffer[i*3 + 2] = 0.8 * (element.cell.energy);
		} else {
			frameBuffer[i*3] = 0;
			frameBuffer[i*3 + 1] = element.environment.chem * 0.6;
			frameBuffer[i*3 + 2] = 0;
		}
	}
}

void VideoOutput::getChemFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx) {
	for (int i = 0; i < nPixels; i++) {
		GridElement element = grid[i + dIdx*nPixels];
		if (element.cell.alive) {
			frameBuffer[i*3] = 0.5 * (255 - element.cell.chem);
			frameBuffer[i*3 + 1] = 0;
			frameBuffer[i*3 + 2] = 0.8 * (element.cell.chem);
		} else {
			frameBuffer[i*3] = 0;
			frameBuffer[i*3 + 1] = element.environment.chem * 0.6;
			frameBuffer[i*3 + 2] = 0;
		}
	}
}

void VideoOutput::getToxinFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx) {
	for (int i = 0; i < nPixels; i++) {
		GridElement element = grid[i + dIdx*nPixels];
		if (element.cell.alive) {
			frameBuffer[i*3] = 0.4 * (255 - element.cell.dToxin);
			frameBuffer[i*3 + 1] = 0.4 * (element.cell.dToxin);
			frameBuffer[i*3 + 2] = 0.4 * (element.cell.dToxin);
		} else {
			frameBuffer[i*3] = (element.environment.ndToxin + element.environment.dToxin) * 0.2;
			frameBuffer[i*3 + 1] = (element.environment.ndToxin + element.environment.dToxin) * 0.2;
			frameBuffer[i*3 + 2] = 0;
		}
	}
}

void VideoOutput::getGenomeFrame(GridElement *grid, int nPixels, unsigned char *frameBuffer, int dIdx) {
	for (int i = 0; i < nPixels; i++) {
		GridElement element = grid[i + dIdx*nPixels];
		if (element.cell.alive && element.cell.genome == 1) {
			frameBuffer[i*3] = 255;
			frameBuffer[i*3 + 1] = 0;
			frameBuffer[i*3 + 2] = 0;
		// } else {
		// 	frameBuffer[i*3] = (element.environment.ndToxin + element.environment.genome) * 0.2;
		// 	frameBuffer[i*3 + 1] = (element.environment.ndToxin + element.environment.dToxin) * 0.2;
		// 	frameBuffer[i*3 + 2] = 0;
		}
		else if(element.cell.alive && element.cell.genome == 2) {
			frameBuffer[i*3] = 0;
			frameBuffer[i*3 + 1] = 0;
			frameBuffer[i*3 + 2] = 255;
		}
		// else {
		// 	frameBuffer[i*3] = 0;
		// 	frameBuffer[i*3 + 1] = 0;
		// 	frameBuffer[i*3 + 2] = 0;
		// }
		else {
			frameBuffer[i*3] = 0;
			frameBuffer[i*3 + 1] = element.environment.chem * 0.6;
			frameBuffer[i*3 + 2] = 0;
		}
	}
}
