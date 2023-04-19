cmake_minimum_required(VERSION "3.16.3")
set (CMAKE_CXX_STANDARD 20)

find_package(Git REQUIRED)

if(NOT EXISTS "${CMAKE_SOURCE_DIR}/external/libsais/LICENSE")
  execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

project (MOVE-r LANGUAGES C CXX)

cmake_policy(SET CMP0081 NEW)
cmake_policy(SET CMP0015 NEW)
cmake_policy(SET CMP0048 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS_DEBUG           "-march=native -DDEBUG -ggdb3 -O0")
set(CMAKE_CXX_FLAGS_RELEASE         "-march=native -DNDEBUG -Ofast")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "-march=native -DDEBUG -ggdb3 -Ofast")
set(BUILD_TESTING OFF)

# openmp
find_package(OpenMP REQUIRED)

# ips4o
add_subdirectory(external/ips4o)
set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

# concurrentqueue
add_subdirectory(external/concurrentqueue/)

# malloc_count
add_library(malloc_count STATIC external/malloc_count/malloc_count.c)
target_include_directories(malloc_count PUBLIC external/malloc_count/)
target_link_libraries(malloc_count PRIVATE dl)

# libsais
add_library(libsais STATIC external/libsais/libsais.c external/libsais/libsais64.c)
target_include_directories(libsais PUBLIC external/libsais/)
target_link_libraries(libsais PRIVATE OpenMP::OpenMP_C)

# sdsl
add_subdirectory(external/sdsl-lite/)
target_compile_options(sdsl PRIVATE -w)
set_target_properties(uninstall-sdsl uninstall gtest gtest_main divsufsort divsufsort64 PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
target_include_directories(sdsl PUBLIC ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include/)

# MOVE-r
add_executable(MOVE-r-genpatterns MOVE-r-genpatterns.cpp)

link_libraries(OpenMP::OpenMP_CXX sdsl malloc_count)
include_directories(internal/)

add_executable(MOVE-r-build MOVE-r-build.cpp)
add_executable(MOVE-r-revert MOVE-r-revert.cpp)
add_executable(MOVE-r-count MOVE-r-count.cpp)
add_executable(MOVE-r-locate MOVE-r-locate.cpp)

target_link_libraries(MOVE-r-build PRIVATE libsais ips4o)
target_include_directories(MOVE-r-build PRIVATE mds/mds_builder/ external/concurrentqueue/)
target_link_libraries(MOVE-r-locate PRIVATE ips4o)