#!/bin/bash

# If Phoenix modules not loaded, load:
modules="CUDA/9.2.148.1 CMake/3.12.1-GCCcore-5.4.0 Anaconda3/2019.03"
for m in $modules; do
	if ! module list 2>&1 | grep -q $m; then
		module load $m
	fi
done
 
# Check if cmake has been run:
if [ ! -d build ] ; then
	mkdir build
	cd build
	cmake ..
	cd ..
fi

# Build:
make -C build -s


# For running:

# Exit upon compilation error:
# if [[ $? -ne 0 ]] ; then
# 	exit 1
# fi

# Run:
# if [[ $# -ne 2 ]] ; then
# 	# Run with default example configuration:
# 	sbatch --wait experiment/job/test_job.sh "experiment/conf/example.json"
# else
# 	sbatch --wait experiment/job/test_job.sh $1
# fi

# Print results:
# cat slurm-*.out