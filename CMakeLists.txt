cmake_minimum_required(VERSION 2.8)
project(treecode_rte)

set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(CMakeScripts_MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_MODULE_PATH ${CMakeScripts_MODULE_DIR})


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -DMATLAB_MEX_FILE -DRUN_OMP -DDISP\
    -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread -DMX_COMPAT_32 \
    -fomit-frame-pointer -malign-double -fstrict-aliasing -march=native -std=gnu++11 -fopenmp -ffast-math")


set(SOURCE_FILES
        blas_wrapper.cpp blas_wrapper.h linalg.h utils.h utility/config.h utility/config.cpp treecode.h treecode_rte.cpp treecode_rte.h main.cpp)

add_executable(TreeCode ${SOURCE_FILES})
target_link_libraries(TreeCode m blas quadmath)