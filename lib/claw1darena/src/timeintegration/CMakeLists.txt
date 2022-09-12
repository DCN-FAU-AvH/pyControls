file(GLOB sources "${PROJECT_SOURCE_DIR}/src/timeintegration/*.cpp")
add_library(timeintegration ${sources})
target_link_libraries(timeintegration grid)
