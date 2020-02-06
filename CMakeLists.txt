cmake_minimum_required(VERSION 3.9)

project(cell_model LANGUAGES CXX CUDA)

add_library(model STATIC
	src/cell_model.cu
	src/cell_model.cuh
	src/cell/movement.cuh
	src/cell/movement.cu
	src/cell_model_params.h
)

target_compile_features(model PUBLIC cxx_std_11)
set_target_properties(model PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(simulate src/main.cpp)
set_property(TARGET simulate PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(simulate PRIVATE model)
