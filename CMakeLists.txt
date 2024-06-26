cmake_minimum_required(VERSION 3.5)

project(cpp_lane_detector LANGUAGES CXX)

option(WITH_VISUALIZATION "Use GUI visualization of lines and point clouds" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(CXX_FILESYSTEM_LIBRARIES stdc++fs)
else()
  set(CXX_FILESYSTEM_LIBRARIES)
endif()

find_package(Eigen3 REQUIRED)

add_executable(cpp_lane_detector
  main.cpp
  point_extractor.cpp
  kdtree.cpp
  clustering.cpp
  lines_detector.cpp
  main_direction_getter.cpp
)
target_link_libraries(cpp_lane_detector PRIVATE Eigen3::Eigen ${CXX_FILESYSTEM_LIBRARIES})

# This part suppresses duplicated branches warning in Eigen library.
if(CMAKE_CXX_COMPILER_ID STREQUAL GNU
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0
)
  target_compile_options(
    cpp_lane_detector
    PUBLIC "-Wno-duplicated-branches"
  )
endif()

if (WITH_VISUALIZATION)
  find_package(PCL 1.8 REQUIRED COMPONENTS common visualization)

  add_library(visualizer_lib OBJECT visualizer.cpp)
  target_include_directories(visualizer_lib PRIVATE ${PCL_INCLUDE_DIRS})
  target_link_libraries(visualizer_lib PUBLIC ${PCL_LIBRARIES})
else()
  add_library(visualizer_lib OBJECT visualizer_stub.cpp)
endif()

target_link_libraries(cpp_lane_detector PRIVATE visualizer_lib)

