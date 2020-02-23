cmake_minimum_required(VERSION 3.9)

project(cell_model LANGUAGES CXX CUDA)

include(FetchContent)

FetchContent_Declare(
	nlohmann_json
	URL https://github.com/nlohmann/json/releases/download/v3.7.3/include.zip
)

FetchContent_Populate(nlohmann_json)

add_library(json INTERFACE)
target_include_directories(json INTERFACE ${nlohmann_json_SOURCE_DIR}/single_include)

# Model:
# TODO: Consider using subdirectory CMakeLists.txt files:
add_library(model STATIC
	src/model/cell_model.cu
	src/model/cell_model.cuh
	src/model/cell_model_params.h
	src/model/cell_model_params.cpp

	# Operations:
	src/model/operation/movement.cuh
	src/model/operation/movement.cu

	# Output:
	src/model/output/video.h
	src/model/output/video.cpp
	src/model/output/statistics.h
	src/model/output/statistics.cpp
)

target_compile_features(model PUBLIC cxx_std_11)
set_target_properties(model PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(model PRIVATE json)

# Driver:
add_executable(simulate src/main.cpp)
set_property(TARGET simulate PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(simulate PRIVATE model)
