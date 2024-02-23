cmake_minimum_required(VERSION 3.5)

project(cpp_lane_detector LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  set(CXX_FILESYSTEM_LIBRARIES stdc++fs)
else()
  set(CXX_FILESYSTEM_LIBRARIES)
endif()

find_package(PCL 1.8 REQUIRED COMPONENTS common io visualization)

add_executable(cpp_lane_detector main.cpp)
target_include_directories(cpp_lane_detector PRIVATE ${PCL_INCLUDE_DIRS})
target_link_libraries(cpp_lane_detector PRIVATE ${PCL_LIBRARIES} ${CXX_FILESYSTEM_LIBRARIES})
