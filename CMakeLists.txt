cmake_minimum_required(VERSION 2.8)

project(smallpt-mpi)

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})
set(LIBRARY_FILES ${LIBRARY_FILES} ${MPI_LIBRARIES})

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})

add_executable(smallpt-mpi src/main.cpp)
target_link_libraries(smallpt-mpi ${LIBRARY_FILES})

