cmake_minimum_required(VERSION 3.9)

project(cell_model LANGUAGES CUDA CXX)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --generate-line-info -g")

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
	src/model/cell_model_driver.cu
	src/model/param/simulation_parameters.cu
	src/model/param/intervention.cu

	# Operations:
	src/model/operation/util.cu
	src/model/operation/movement.cu
	src/model/operation/growth.cu

	# Output:
	src/model/output/video.cu
	src/model/output/statistics.cu
)


target_compile_features(model PUBLIC cxx_std_11)
set_target_properties(model PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(model PRIVATE json)

# Set nvcc flags:
# target_compile_options(model PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
#    --generate-line-info
#    -g
# >)

set_property(TARGET model PROPERTY CUDA_ARCHITECTURES 61-real 61-virtual)


# Driver:
add_executable(simulate src/main.cpp)
set_property(TARGET simulate PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(simulate PRIVATE model)	
