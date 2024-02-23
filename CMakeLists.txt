cmake_minimum_required(VERSION 3.5)

project(cpp_lane_detector LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  set(CXX_FILESYSTEM_LIBRARIES stdc++fs)
else()
  set(CXX_FILESYSTEM_LIBRARIES)
endif()

add_executable(cpp_lane_detector main.cpp point_extractor.cpp)
target_link_libraries(cpp_lane_detector PRIVATE ${CXX_FILESYSTEM_LIBRARIES})


find_package(PCL 1.8 REQUIRED COMPONENTS common visualization)

add_library(visualizer_lib STATIC visualizer.cpp)
target_include_directories(visualizer_lib PRIVATE ${PCL_INCLUDE_DIRS})
target_link_libraries(visualizer_lib PUBLIC ${PCL_LIBRARIES})

target_link_libraries(cpp_lane_detector PRIVATE visualizer_lib)

