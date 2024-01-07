cmake_minimum_required( VERSION 3.16 )
project(disjoint_path)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -g -Wall -rdynamic -pthread -march=native")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

set(MAIN_FILE_LIST
    tools.cpp
    Graph.cpp
    Baseline.cpp
    BatchSet.cpp
    MultiVP.cpp

    main.cpp
)
add_executable(vdp ${MAIN_FILE_LIST})

set(DATASET_FILE_LIST
    tools.cpp
    Graph.cpp

    dataset.cpp
)
add_executable(dataset ${DATASET_FILE_LIST})