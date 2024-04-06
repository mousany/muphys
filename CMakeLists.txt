cmake_minimum_required(VERSION 3.20)
project(muphys_cpp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
# notify all possible issues 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

option(MU_ENABLE_SINGLE "Enable single precision" OFF)
# available front-ends
option(MU_IMPL "Select implementation" "seq")

option(MU_ENABLE_TESTS "Enable unit-tests" ON)

if(NOT MU_GPU_VERSION)
    set(MU_GPU_VERSION "cc70")
endif()

# includes
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# configure the structure of the output folders
include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}" CACHE PATH
        "Directory for the built binaries")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
        "Directory for the built libraries")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
        "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}" CACHE PATH
        "Directory for the built static libraries")

if (MU_ENABLE_SINGLE)
    add_compile_definitions(__SINGLE_PRECISION)
endif ()

if(NOT MU_IMPL)
    message(FATAL_ERROR "No implementation was selected. Options are: seq, omp, gpu")
# add compile definitions and external dependecies based on options
elseif ("${MU_IMPL}" STREQUAL "seq")
    add_compile_definitions(MU_ENABLE_SEQ)
elseif("${MU_IMPL}" STREQUAL "omp")
    add_compile_definitions(MU_ENABLE_OMP)
    add_compile_options(-march=native)
elseif ("${MU_IMPL}" STREQUAL "gpu")
    add_compile_definitions(MU_ENABLE_GPU)
    #set(CMAKE_C_COMPILER nvc)
    #set(CMAKE_CXX_COMPILER nvc++)
    if (NOT CLION)
        add_compile_options(-mp=gpu -gpu=${MU_GPU_VERSION})
        add_link_options(-mp=gpu -gpu=${MU_GPU_VERSION})
        #target_sources(muphys_cpp PRIVATE "implementations/gpu/graupel.cpp")
    endif()
else()
    message(FATAL_ERROR "${MU_IMPL} is not a valid configuration.")
endif ()

# add local sources
add_subdirectory(core)
add_subdirectory(io)
include(implementations/CMakeLists.txt)

# add test (if enabled)
if (MU_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(extern/googletest)
    add_subdirectory(test)
endif ()

# build project executable
add_executable(graupel "main.cpp")

# link dependency libs
target_link_libraries(graupel muphys_core muphys_io muphys_implementation)
target_include_directories(graupel PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/core"
                          "${PROJECT_SOURCE_DIR}/implementations"
                          "${PROJECT_SOURCE_DIR}/io"
                        )
set_target_properties(graupel PROPERTIES LINKER_LANGUAGE CXX)
