cmake_minimum_required(VERSION 3.3)
project(hpxct CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HPX_IGNORE_COMPILER_COMPATIBILITY ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#add_definitions(-DENABLE_DEBUG_LOGGING)
#add_definitions(-DHPXIC_ENABLE_APEX=ON)
add_definitions(-DFLATAUGMENTATION)
#add_definitions(-DFILEOUT)
add_definitions(-DVTIOUT)
#add_definitions(-DLAZYREMOTECALLS)
#add_definitions(-DENABLE_APEX_PROFILING)

find_package(HPX REQUIRED)
find_package(VTK REQUIRED)
find_package(Boost REQUIRED)
find_package(glm REQUIRED)

include(${VTK_USE_FILE})

include_directories(${HPX_INCLUDE_DIR} ${GLM_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

add_hpx_component(TreeConstructor
    ESSENTIAL
    SOURCES TreeConstructor.cpp DataManager.cpp
)

add_hpx_executable(hpxct
    ESSENTIAL
    SOURCES main.cpp 
    COMPONENT_DEPENDENCIES TreeConstructor
    DEPENDENCIES ${Boost_LIBRARIES} ${VTK_LIBRARIES} ${TEEM_LIBRARIES}
)


