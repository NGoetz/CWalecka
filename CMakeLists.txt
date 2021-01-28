cmake_minimum_required(VERSION 3.16)
project(cwalecka0.1
    VERSION 0.1
    DESCRIPTION "C++ Walecka model implementation"
    LANGUAGES CXX)
add_subdirectory(matplotplusplus)
add_executable(walecka main.cpp
    src/wplot.cpp
    src/wplot.h 
    src/wsolvers.cpp 
    src/wsolvers.h 
    src/wfunctions.cpp 
    src/wfunctions.h )
find_package(GSL REQUIRED)    
target_link_libraries(walecka PUBLIC matplot GSL::gsl GSL::gslcblas)
target_compile_features(walecka PRIVATE cxx_std_20)
