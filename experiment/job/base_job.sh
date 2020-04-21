#!/bin/bash

#SBATCH -p gpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=00:20:00
#SBATCH --gres=gpu:kepler:1
#SBATCH --mem=1000MB

module load CUDA/9.2.148.1
module load FFmpeg/3.3.1-foss-2016b

./build/simulate $1
