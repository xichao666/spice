cmake_minimum_required(VERSION 3.20)

project(spice_dc_solver LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(SPICE_BUILD_LEGACY_EXAMPLES
    "Build legacy hand-written BJT examples that have not yet migrated to DcProblem" OFF)

set(DC_SOLVER_SOURCES
    src/dc_solver.c
    src/homotopy.c
    src/gmin_stepping.c
    src/source_stepping.c
    src/transient_solver.c
)

add_library(dcsolve STATIC ${DC_SOLVER_SOURCES})

target_include_directories(dcsolve PUBLIC include)

if(MSVC)
    target_compile_options(dcsolve PRIVATE /W4)
else()
    target_compile_options(dcsolve PRIVATE -Wall -Wextra -Wpedantic)
    target_link_libraries(dcsolve PUBLIC m)
endif()

add_executable(netlist_dc_solver src/netlist_dc_solver.c)
target_link_libraries(netlist_dc_solver PRIVATE dcsolve)

if(SPICE_BUILD_LEGACY_EXAMPLES)
    foreach(circuit IN ITEMS bjt1 bjt2 bjt3)
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${circuit}_source_stepping.c")
            add_executable(${circuit}_source_stepping
                examples/${circuit}_source_stepping.c)
            target_link_libraries(${circuit}_source_stepping PRIVATE dcsolve)
        endif()

        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${circuit}_nr.c")
            add_executable(${circuit}_nr examples/${circuit}_nr.c)
            target_link_libraries(${circuit}_nr PRIVATE dcsolve)
        endif()
    endforeach()
endif()

enable_testing()

add_executable(test_dc_solver tests/test_dc_solver.c)
target_link_libraries(test_dc_solver PRIVATE dcsolve)
add_test(NAME dc_solver_smoke COMMAND test_dc_solver)

if(SPICE_BUILD_LEGACY_EXAMPLES AND
   EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/bjt1_source_stepping.c")
    add_executable(source_stepping_benchmark tests/source_stepping_benchmark.c)
    target_link_libraries(source_stepping_benchmark PRIVATE dcsolve)
endif()

add_executable(test_homotopy tests/test_homotopy.c)
target_link_libraries(test_homotopy PRIVATE dcsolve)
add_test(NAME fixed_point_homotopy_smoke COMMAND test_homotopy)

add_executable(reference_homotopy_experiments tests/reference_homotopy_experiments.c)
target_link_libraries(reference_homotopy_experiments PRIVATE dcsolve)
add_test(NAME reference_fixed_point_homotopy COMMAND reference_homotopy_experiments)

add_executable(test_transient_gear tests/test_transient_gear.c)
target_link_libraries(test_transient_gear PRIVATE dcsolve)
add_test(NAME transient_gear COMMAND test_transient_gear)
