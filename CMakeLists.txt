cmake_minimum_required(VERSION 3.16)
project(cwalecka0.1
    VERSION 0.1
    DESCRIPTION "C++ Walecka model implementation"
    LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 20)
list(APPEND flags  "-Wall")
add_subdirectory(thirdparty/matplotplusplus)
set(SOURCES main.cpp src/wplot.h src/wplot.cpp src/wfunctions.cpp src/wfunctions.h src/wsolvers.cpp src/wsolvers.h  )
add_executable(walecka ${SOURCES})
target_compile_options(walecka PRIVATE ${flags})
target_include_directories(walecka PRIVATE src)
target_include_directories(walecka PRIVATE thirdparty/matplotplusplus/source)
find_package(GSL REQUIRED)  
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
target_link_libraries(walecka  PUBLIC Eigen3::Eigen PUBLIC matplot GSL::gsl GSL::gslcblas PRIVATE OpenMP::OpenMP_CXX)


target_compile_features(walecka PRIVATE cxx_std_20)
