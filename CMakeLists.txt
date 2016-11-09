project(MapBuilder)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH "/usr/local/lib/cmake/" )
#find_package(G2O REQUIRED)
find_package(YARP REQUIRED)
FIND_PACKAGE(ICUB REQUIRED)
find_package( OpenCV REQUIRED )
find_package(PCL 1.7 REQUIRED)
FIND_PACKAGE(SuiteSparse REQUIRED)
find_package(Ceres REQUIRED)
find_library(CHOLMOD_LIBRARY NAMES cholmod)
list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4")
# add YARP include directories
include_directories(${YARP_INCLUDE_DIRS}
                    ${ICUB_INCLUDE_DIRS}
                    ${OpenCV_INCLUDE_DIRS}
                    ${PCL_INCLUDE_DIRS}
                    #${G2O_INCLUDE_DIR}
                    ${CHOLMOD_INCLUDE_DIR}
                    ${SUITESPARSE_INCLUDE_DIRS}
                    ${CERES_INCLUDE_DIRS})
aux_source_directory(. SRC_LIST)

set(FAST_SRCS fast_10.c
              fast_11.c
              fast_12.c
              fast_9.c
              fast.c
              nonmax.c)

set(FAST_HDRS fast.h)

add_library(FAST STATIC ${FAST_SRCS} ${FAST_HDRS})
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(MapBuilder_SRCS featurefinderandtracker.cpp
                    icubimgkingrabber.cpp
                    main.cpp
                    MatrixVecAlgebra.cpp
                    triangulationclass.cpp
                    visibilitymatrix.cpp)

set(MapBuilder_HDRS bal_problem.h
                    featurefinderandtracker.h
                    icubimgkingrabber.h
                    MatrixVecAlgebra.h
                    my_reprojection_error.h
                    triangulationclass.h
                    visibilitymatrix.h)

add_executable(${PROJECT_NAME} ${MapBuilder_SRCS} ${MapBuilder_HDRS})
target_link_libraries(MapBuilder ${YARP_LIBRARIES}
                                         ${ICUB_LIBRARIES}
                                         ${OpenCV_LIBRARIES}
                                         ${PCL_LIBRARIES}
                                         #${G2O_LIBRARIES}
                                         #${G2O_CORE_LIBRARY}
                                         #${G2O_STUFF_LIBRARY}
                                         #${G2O_SOLVER_CHOLMOD}
                                         ${CHOLMOD_LIBRARY}
                                         #${G2O_TYPES_SBA}
                                         #${G2O_TYPES_ICP}
                                        ${CERES_LIBRARIES}
                                    FAST)


if(${CMAKE_VERSION} VERSION_LESS 3.1)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
  check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  endif()
else()
  target_compile_features(MapBuilder PRIVATE cxx_range_for)
endif()
