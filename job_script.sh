#!/bin/bash

#SBATCH -p test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=00:02:00
#SBATCH --gres=gpu:1
#SBATCH --mem=50MB

module load Anaconda3/5.0.1
source activate cell

python main.py

source deactivate