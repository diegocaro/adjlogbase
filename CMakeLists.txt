########################################################################
# CMake build script for libcds.


########################################################################
#
# Project-wide settings

# Name of the project.
#
# CMake files in this project can refer to the root source directory
# as ${cds_SOURCE_DIR} and to the root binary directory as
# ${cds_BINARY_DIR}.


cmake_minimum_required(VERSION 2.8.7)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(AppendCompilerFlags)

project(evelog CXX)

### Flags for all
SET(CMAKE_CXX_FLAGS "-Wall -DEXPERIMENTS")

# C++11 compiler Check
if(NOT CMAKE_CXX_COMPILER_VERSION) # work around for cmake versions smaller than 2.8.10
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
endif()
if(CMAKE_CXX_COMPILER MATCHES ".*clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(CMAKE_COMPILER_IS_CLANGXX 1)
endif()
if( (CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.7) OR
    (CMAKE_COMPILER_IS_CLANGXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.2))
  message(FATAL_ERROR "Your C++ compiler does not support C++11. Please install g++ 4.7 (or greater) or clang 3.2 (or greater)")
else()
  message(STATUS "Compiler is recent enough to support C++11.")
endif()

if( CMAKE_COMPILER_IS_GNUCXX )
    append_cxx_compiler_flags("-std=c++11 -Wall -Wextra  -DNDEBUG" "GCC" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_OPT_FLAGS)
    if ( CODE_COVERAGE )
        append_cxx_compiler_flags("-g -fprofile-arcs -ftest-coverage -lgcov" "GCC" CMAKE_CXX_FLAGS)
    endif()
else()
    append_cxx_compiler_flags("-std=c++11 -DNDEBUG" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "CLANG" CMAKE_CXX_OPT_FLAGS)
endif()




## lib for compressing small integers
add_subdirectory(../compresslists/ ext_compresslists)

add_subdirectory(../etdc/ ext_etdc)

include_directories(
  ${evelog_SOURCE_DIR}
)

add_library(evelog tgraph.cpp arraysort.cpp)
target_link_libraries(evelog compresslists etdc)

add_executable(create create.cpp )
target_link_libraries(create evelog)

add_executable(use use.cpp )
target_link_libraries(use evelog)

add_executable(benchmark benchmark.cpp )
target_link_libraries(benchmark evelog)

add_executable(benchmark2 benchmark2.cpp )
target_link_libraries(benchmark2 evelog)


add_executable(exp-query exp-query.cpp )
target_link_libraries(exp-query evelog)

add_executable(vertices vertices.cpp )
target_link_libraries(vertices evelog)



if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(benchmark -lrt)
target_link_libraries(benchmark2 -lrt)
target_link_libraries(exp-query -lrt)
target_link_libraries(vertices -lrt)
endif ()