#!/bin/bash

if [ ! -d build ] ; then
	mkdir build
	cd build
	cmake ..
	cd ..
fi

# build:
rm slurm-*.out
make -C build -s
# Check for compilation error:
if [[ $? -ne 0 ]] ; then
	exit 1
fi

# run:
sbatch test_job.sh

# output results:
sleep 3
cat slurm-*.out
sleep 9
cat slurm-*.out