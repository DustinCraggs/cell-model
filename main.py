import argparse
import numpy as np
import numba
from numba import cuda
from numba.cuda import random

from util.save_video import video_stream

# Simulation parameters:
x_dim, y_dim = 250, 250
init_prob = 0.1
energy_survival_threshold = 10

light_function = lambda _pos: 3
gather_light_energy_cost = 1
gather_light_energy_prob = 0.3

co2_function = lambda _pos: 1

# Data packing: (TODO: calc offset)
energy_offset = 1
energy_bits = 6
chem_offset = 7
chem_bits = 6

# Cuda parameters:
# Threads per block:
tpb = 4
n_threads = (tpb, tpb)
n_blocks = ((x_dim + tpb - 1) // tpb, (y_dim + tpb - 1) // tpb)

def main(args):
	cells = allocate_cell_memory()
	rng_states = allocate_rng_states()
	initialise_cells[n_blocks, n_threads](cells, rng_states, init_prob)
	with video_stream((x_dim, y_dim), 'occu.mp4', scale=3) as occu_out, \
		 video_stream((x_dim, y_dim), 'enrg.mp4', scale=3) as enrg_out:
		for i in range(100):
			print(i, end=', ', flush=True)
			simulate_cells[n_blocks, n_threads](cells, rng_states, 1)
			occu_out.write(occupied(np.array(cells)) * 255)
			enrg_out.write(energy(np.array(cells)) * 4)
		print()

def print_state(cells):
	print("Occupation:")
	print(occupied(np.array(cells)))
	print("Energy:")
	print(energy(np.array(cells)))

def allocate_cell_memory():
	return cuda.device_array((x_dim, y_dim), dtype=np.int32)

def allocate_rng_states():
	return random.create_xoroshiro128p_states(x_dim * y_dim, seed=1)

@cuda.jit
def initialise_cells(cells, rng_states, prob):
	i = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
	j = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
	if i < cells.shape[0] and j < cells.shape[1]:
		cells[i][j] = 0
		rng_i = i + j * x_dim
		if random.xoroshiro128p_uniform_float32(rng_states, rng_i) < prob:
			cells[i][j] = 1
			energy = int(random.xoroshiro128p_uniform_float32(rng_states, rng_i) * 2**energy_bits)
			cells[i][j] += energy << energy_offset

@cuda.jit
def simulate_cells(cells, rng_states, iterations):
	i = cuda.blockIdx.x * cuda.blockDim.x + cuda.threadIdx.x
	j = cuda.blockIdx.y * cuda.blockDim.y + cuda.threadIdx.y
	if i < cells.shape[0] and j < cells.shape[1]:
		for _ in range(iterations):
			iterate_cells(cells, rng_states, i, j)

@cuda.jit(device=True, inline=True)
def iterate_cells(cells, rng_states, i, j):
	# All kernels/ops for one iteration here:
	if occupied(cells[i][j]):
		check_survival(cells, i, j)
		reduce_energy(cells, i, j)
		acquire_energy(cells, rng_states, i, j)

@cuda.jit(device=True, inline=True)
def reduce_energy(cells, i, j):
	if energy(cells[i][j]) > 0:
		set_energy(cells, i, j, energy(cells[i][j]) - 1)

@cuda.jit(device=True, inline=True)
def acquire_energy(cells, rng_states, i, j):
	if 0 < energy(cells[i][j]) < 2**energy_bits:
		rng_i = i + j * x_dim
		if random.xoroshiro128p_uniform_float32(rng_states, rng_i) < gather_light_energy_prob:
			new_energy = energy(cells[i][j]) + light_function((i, j)) - gather_light_energy_cost
			new_energy = min(new_energy, 2**energy_bits - 1)
			set_energy(cells, i, j, new_energy)

@cuda.jit(device=True, inline=True)
def check_survival(cells, i, j):
	if energy(cells[i][j]) < energy_survival_threshold:
		cells[i][j] = 0

@cuda.jit(device=True, inline=True)
def occupied(cell):
	return cell & 1

@cuda.jit(device=True, inline=True)
def energy(cell):
	return cell >> energy_offset & (0xffffffff >> (32 - energy_bits))

@cuda.jit(device=True, inline=True)
def set_energy(cells, i, j, val):
	inv_energy_mask = ((0xffffffff >> (32 - energy_offset)) | (0xffffffff << (energy_bits + energy_offset)))
	cells[i][j] = cells[i][j] & inv_energy_mask
	cells[i][j] += val << energy_offset

def chem(cell):
	return cell >> chem_offset & (0xffffffff >> (32 - chem_bits))

def _parse_args():
	parser = argparse.ArgumentParser(description='Cell test')
	return parser.parse_args()

if __name__ == '__main__':
	main(_parse_args())