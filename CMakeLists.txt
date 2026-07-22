cmake_minimum_required(VERSION 3.20)

project(spice_dc_solver LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(dcsolve STATIC
    src/dc_solver.c
    src/source_stepping.c
)

target_include_directories(dcsolve PUBLIC include)

if(MSVC)
    target_compile_options(dcsolve PRIVATE /W4)
else()
    target_compile_options(dcsolve PRIVATE -Wall -Wextra -Wpedantic)
    target_link_libraries(dcsolve PUBLIC m)
endif()

add_executable(netlist_dc_solver src/netlist_dc_solver.c)
target_link_libraries(netlist_dc_solver PRIVATE dcsolve)

foreach(circuit IN ITEMS bjt1 bjt2 bjt3)
    add_executable(${circuit}_source_stepping
        examples/${circuit}_source_stepping.c)
    target_link_libraries(${circuit}_source_stepping PRIVATE dcsolve)

    add_executable(${circuit}_nr examples/${circuit}_nr.c)
    target_link_libraries(${circuit}_nr PRIVATE dcsolve)
endforeach()

enable_testing()

add_executable(test_dc_solver tests/test_dc_solver.c)
target_link_libraries(test_dc_solver PRIVATE dcsolve)
add_test(NAME dc_solver_smoke COMMAND test_dc_solver)
