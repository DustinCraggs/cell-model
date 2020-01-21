#!/bin/bash

# build:
rm slurm-*.out
make -C build -s

# run:
sbatch test_job.sh

# output results:
sleep 3
cat slurm-*.out
sleep 9
cat slurm-*.out