cmake_minimum_required(VERSION 3.22)
PROJECT(DemoDGtal)

FIND_PACKAGE(DGtal REQUIRED)

ADD_EXECUTABLE(main main.cpp)
TARGET_LINK_LIBRARIES(main DGtal::DGtal)
