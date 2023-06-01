cmake_minimum_required(VERSION 3.1)
project(derna)

set(CMAKE_CXX_STANDARD 11)

IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE Release)
ENDIF()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS "-Wno-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


add_executable(derna src/main.cpp src/utils.h src/utils.cpp src/Nussinov.cpp src/Nussinov.h src/NussinovAlgorithm.cpp src/NussinovAlgorithm.h src/Zuker.cpp src/Zuker.h src/default.cpp src/default.h src/params/constants.h src/params/intl11.h src/params/intl11dH.h src/params/intl21.h src/params/intl21dH.h src/params/intl22.h src/params/intl22dH.h src/ZukerAlgorithm.cpp src/ZukerAlgorithm.h)

