# All test executables: link to flux_physics (defined in parent)
add_executable(testFlux testFlux.cpp)
target_link_libraries(testFlux PRIVATE flux_physics)
add_custom_target(run_test COMMAND testFlux DEPENDS testFlux)

add_executable(test_mesh test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE mesh_io)
add_custom_target(run_test_mesh COMMAND test_mesh DEPENDS test_mesh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_verify test_verify.cpp)
target_link_libraries(test_verify PRIVATE mesh_io bc)
add_custom_target(run_test_verify COMMAND test_verify DEPENDS test_verify
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(test_bc_flux test_bc_flux.cpp)
target_link_libraries(test_bc_flux PRIVATE flux_physics)

# Add more tests the same way:
# add_executable(testXxx testXxx.cpp)
# target_link_libraries(testXxx PRIVATE flux_physics)
# add_custom_target(run_test_xxx COMMAND testXxx DEPENDS testXxx)
