cmake_minimum_required(VERSION 3.9)
project(derna)
cmake_policy(SET CMP0069 NEW)

set(CMAKE_CXX_STANDARD 17)

IF( NOT CMAKE_BUILD_TYPE )
    SET( CMAKE_BUILD_TYPE Release)
ENDIF()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-uninitialized")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# OpenMP support (needed for Block 1 parallelization in PositionBeamDP).
# On macOS with Apple Clang, libomp from Homebrew is required.
find_package(OpenMP COMPONENTS CXX QUIET)
if(NOT OpenMP_CXX_FOUND AND APPLE)
    execute_process(COMMAND brew --prefix libomp
                    OUTPUT_VARIABLE LIBOMP_PREFIX
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    if(LIBOMP_PREFIX)
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
        find_package(OpenMP COMPONENTS CXX REQUIRED)
    endif()
endif()

add_executable(derna src/main.cpp src/utils.h src/utils.cpp src/Nussinov.cpp src/Nussinov.h src/NussinovAlgorithm.cpp src/NussinovAlgorithm.h src/Zuker.cpp src/Zuker.h src/PositionBasedBeamZuker.cpp src/PositionBasedBeamZuker.h src/PositionBeamDP.cpp src/PositionBeamDP.h src/default.cpp src/default.h src/params/constants.h src/params/intl11.h src/params/intl11dH.h src/params/intl21.h src/params/intl21dH.h src/params/intl22.h src/params/intl22dH.h src/ZukerAlgorithm.cpp src/ZukerAlgorithm.h)
if(OpenMP_CXX_FOUND)
    target_link_libraries(derna PUBLIC OpenMP::OpenMP_CXX)
endif()
