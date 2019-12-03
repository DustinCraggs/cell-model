#!/bin/bash

#SBATCH -p test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=00:01:00
#SBATCH --gres=gpu:1
#SBATCH --mem=5MB

module load CUDA/9.2.148.1

./main