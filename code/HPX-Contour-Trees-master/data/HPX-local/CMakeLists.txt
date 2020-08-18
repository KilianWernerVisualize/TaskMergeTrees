cmake_minimum_required(VERSION 3.3)
project(hpxct-local CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HPX_IGNORE_COMPILER_COMPATIBILITY ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(HPX REQUIRED)
find_package(VTK REQUIRED)
find_package(Boost REQUIRED)
find_package(Teem REQUIRED)

include(${VTK_USE_FILE})

include_directories(libs/ConcurrentQueue ${HPX_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${TEEM_INCLUDE_DIR})


#add_hpx_component(BlockContainer
#    ESSENTIAL
#    SOURCES Container.cpp Block.cpp BlockMap.cpp Image.cpp ImagePool.cpp Compression.cpp Renderer.cpp DataReader.cpp BlockCache.cpp BlockCache.h Container.h Block.h BlockMap.h BlockPriorityScheduler.cpp BlockPriorityScheduler.h Log.h Image.h ImagePool.h Compression.h Serialization.h Config.h Renderer.h DataReader.h NrrdReader.h MandelbulbReader.h CommunicationHandler.h CommunicationHandler.cpp AABB.h Frustum.h Monitor.h
#)

add_hpx_executable(hpxct-local
    ESSENTIAL
    SOURCES main.cpp Arc.h SweepQueue.h Boundary.h DataManager.h Visualizer.h TreeConstructor.h TreeConstructor.cpp Augmentation.h
    #COMPONENT_DEPENDENCIES BlockContainer
    DEPENDENCIES ${Boost_LIBRARIES} ${VTK_LIBRARIES} ${TEEM_LIBRARIES}
)


