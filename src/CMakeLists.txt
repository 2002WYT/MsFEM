cmake_minimum_required(VERSION 3.15)
include_directories(../include)
include_directories(
  ${catkin_INCLUDE_DIRS}
  /opt/intel/oneapi/mkl/2024.1/include/
)
aux_source_directory(. SRC_SOURCE)
add_library(SRC ${SRC_SOURCE})