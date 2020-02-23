#include "video.h"

#include <string>
#include <stdlib.h>

using std::string;

VideoOutput::VideoOutput(CellModelParams params, OutputParams outputParams) {
	energyPipe = openVideoOutputProcess(
		outputParams.video.energyPath,
		params.w,
		params.h
	);
	chemPipe = openVideoOutputProcess(
		outputParams.video.chemPath,
		params.w,
		params.h
	);
	toxinPipe = openVideoOutputProcess(
		outputParams.video.toxinPath,
		params.w,
		params.h
	);

	std::cout << "EP " << energyPipe << std::endl;
	std::cout << "CP " << chemPipe << std::endl;
	std::cout << "TP " << toxinPipe << std::endl;

	nPixels = params.w * params.h;
	frameBuffer = new unsigned char[nPixels * 3];
}

void VideoOutput::writeFrame(CellModel model) {
	GridElement *grid;
	if (energyPipe || chemPipe || toxinPipe) {
		grid = model.getHostGrid();
	}

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
}

void VideoOutput::close() {
	// TODO: Check if opened
	// if (!pclose(energyPipe) || !pclose(chemPipe) || !pclose(toxinPipe)) {
	// 	std::cerr << "Error: Failed to close one or more FFmpeg output processes" << std::endl;
	// }
	if (pclose(energyPipe) != 0) {
		std::cerr << "Error: Failed to close energy FFmpeg output process" << std::endl;
	}
	if (pclose(chemPipe) != 0) {
		std::cerr << "Error: Failed to close chem FFmpeg output process" << std::endl;
	}
	if (pclose(toxinPipe) != 0) {
		std::cerr << "Error: Failed to close toxin FFmpeg output process" << std::endl;
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
		+ path
		+ string(" > debug/ffmpeg.txt");
	std::cout << ffmpegCommand << std::endl;
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
