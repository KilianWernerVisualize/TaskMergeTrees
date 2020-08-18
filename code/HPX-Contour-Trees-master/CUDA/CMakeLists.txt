cmake_minimum_required(VERSION 3.8)
project(cudact CXX CUDA)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CUDA_FLAGS "--gpu-architecture=compute_61")
set(CMAKE_CUDA_LINK_FLAGS "-L /usr/lib -lhalloc")

add_definitions(-DNO_HPX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(VTK REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED )

include(${VTK_USE_FILE})

include_directories(/usr/include ${GLM_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_executable(cudact main.cpp TreeConstructor.cu)

target_compile_features(cudact PUBLIC cxx_std_14)

target_link_libraries(cudact /usr/lib/libhalloc.a ${GLM_LIBRARIES} ${Boost_LIBRARIES} ${VTK_LIBRARIES})




