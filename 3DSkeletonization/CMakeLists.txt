cmake_minimum_required(VERSION 2.8)
 
PROJECT(3DSkeletonization)

find_package(VTK REQUIRED)


INCLUDE_DIRECTORIES(Data)

add_subdirectory(Code)

enable_testing()
