cmake_minimum_required(VERSION 3.14)
project(GjkEpa)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # Export hints for clangd

# Required packages
find_package(Eigen3 REQUIRED)

# Build library
add_library(gjkepa src/gjk.cpp src/epa.cpp)
target_include_directories(gjkepa PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(gjkepa PUBLIC Eigen3::Eigen)

# Build unit tests
option(BUILD_TESTS "Build unit tests" ON)
if(BUILD_TESTS)
    target_compile_definitions(gjkepa PUBLIC BUILD_TESTS)
    add_subdirectory(tests)
endif()
