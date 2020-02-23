#!/bin/bash

#SBATCH -p test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=00:01:00
#SBATCH --gres=gpu:1
#SBATCH --mem=20MB

module load CUDA/9.2.148.1
module load FFmpeg/3.3.1-foss-2016b

mkdir -p experiment/results

./build/simulate $1

# ./build/simulate > results/statistics.csv
# ./build/simulate | ffmpeg -hide_banner -loglevel panic -y -f rawvideo -pixel_format rgb24 -video_size 128x128 \
# 	-i - -c:v h264 -pix_fmt yuv420p -s 512x512 -sws_flags neighbor results/video.mp4

# ffmpeg -hide_banner -loglevel panic -y -f rawvideo -pixel_format rgb24 -video_size 128x128 \
# 	-i - -c:v h264 -pix_fmt yuv420p -s 512x512 -sws_flags neighbor results/video.mp4

