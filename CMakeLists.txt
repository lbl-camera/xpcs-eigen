project(xpcs-eigen)

cmake_minimum_required(VERSION 3.5)

set(CMAKE_BUILD_TYPE RELEASE)
#set(CMAKE_BUILD_TYPE DEBUG)

# build python module
option(ALS_BUILD "build python module instead of executable" OFF)
if (ALS_BUILD) 
    message(STATUS "builing python library instead of executable")
endif ()

# c++ stuff
add_definitions(-std=c++0x)

# sources 
set(EXESRC src/main.cpp)
set(LIBSRC src/py_xpcs.cpp)


add_subdirectory(src/xpcs)
add_subdirectory(src/io)

find_package(OpenMP)
if (OPENMP_FOUND)
    list(APPEND ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
    list(APPEND ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
endif()

# hdf5
find_package(HDF5)
if (HDF5_FOUND) 
    add_definitions(-DHAVE_HDF5)
    set(LIBS ${LIBS} ${HDF5_LIBRARIES})
    include_directories(${HDF5_INCLUDE_DIRS})
endif()

# eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# gflags
find_package(gflags REQUIRED)

# spdlog
find_package(spdlog)

# boost python
if (ALS_BUILD)
    # python library
    find_package(PythonLibs 3.7 REQUIRED)
    list(APPEND ${LIBS} ${PYTHON_LIBRARIES})
    include_directories(${PYTHON_INCLUDE_DIRS})

    # boost.python
    find_package(Boost REQUIRED COMPONENTS python37 numpy37)
    add_definitions(-DHAVE_BOOST_PYTHON)
    list(APPEND ${LIBS} ${Boost_LIBRARIES})
    include_directories(${Boost_INCLUDE_DIRS})
endif (ALS_BUILD)


if (ALS_BUILD)
    message(STATUS ${CMAKE_CXX_FLAGS})
    add_library(pyXPCS SHARED ${LIBSRC})
    target_link_libraries(pyXPCS ${LIBS})
    set_target_properties(pyXPCS PROPERTIES
                                 PREFIX ""
                         )
else (ALS_BUILD)
    include_directories("src" "src/external")
    add_executable(corr ${EXESRC})
    target_link_libraries(corr hdf5 gflags)
endif(ALS_BUILD)
