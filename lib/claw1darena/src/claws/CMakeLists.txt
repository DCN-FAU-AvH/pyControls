set (ENABLE_WARNINGS ${VERBOSE_MAKE})
if (VERBOSE_MAKE)
	set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

add_library(claws INTERFACE)
target_include_directories(claws INTERFACE ${PROJECT_SOURCE_DIR}/src/claws/)
