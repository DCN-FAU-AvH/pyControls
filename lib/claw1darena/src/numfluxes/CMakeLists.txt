file(GLOB sources "${PROJECT_SOURCE_DIR}/src/numfluxes/*.cpp")
add_library(numfluxes ${sources})
target_link_libraries(numfluxes dof)
