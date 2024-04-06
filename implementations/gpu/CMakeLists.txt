add_library(muphys_implementation SHARED "graupel.cpp")
target_link_libraries(muphys_implementation muphys_core)

find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(muphys_implementation OpenMP::OpenMP_CXX)
else()
    message(FATAL_ERROR "OpenMP not found")
endif()
