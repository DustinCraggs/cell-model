#!/bin/bash

# Check if cmake has been run:
if [ ! -d build ] ; then
	mkdir build
	cd build
	cmake ..
	cd ..
fi

# Move any previous results to results/old/
mkdir -p results/old/
mv slurm-*.out results/old/

# Build:
make -C build -s
# Exit upon compilation error:
if [[ $? -ne 0 ]] ; then
	exit 1
fi

# Run:
sbatch --wait test_job.sh

# Print results:
cat slurm-*.out