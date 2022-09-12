file(GLOB sources "${PROJECT_SOURCE_DIR}/src/utils/*.cpp")
add_library(utils ${sources})
target_link_libraries(utils grid)
