# Module that checks whether the compiler supports the
# quadruple precision floating point math
#
# Sets the following variables:
# HAVE_QUADMATH
# QUADMATH_LIBRARIES
#
# perform tests
include(CheckCXXCompilerFlag)
include(CheckCSourceCompiles)
include(CMakePushCheckState)

if(USE_QUADMATH)
  if(NOT DEFINED HAVE_EXTENDED_NUMERIC_LITERALS)
    check_cxx_compiler_flag("-Werror -fext-numeric-literals" HAVE_EXTENDED_NUMERIC_LITERALS)
  endif()

  if (HAVE_EXTENDED_NUMERIC_LITERALS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
  endif()

  cmake_push_check_state(RESET)
  list(APPEND CMAKE_REQUIRED_LIBRARIES "quadmath")
  CHECK_CXX_SOURCE_COMPILES("
#include <quadmath.h>
int main(void){
    __float128 x = sqrtq(1.0);
    x = FLT128_MIN;
}" QUADMATH_FOUND)
  cmake_pop_check_state()

  #(at least) gcc7 already defines the overload for std::abs<__float128>
  cmake_push_check_state(RESET)
  CHECK_CXX_SOURCE_COMPILES("
#include <quadmath.h>
int main(void){
    __float128 x = 1.0;
    x = std::abs(1.0);
}" STD_HAVE_ABSQUAD)
  cmake_pop_check_state()

  if (QUADMATH_FOUND)
    set(QUADMATH_LIBRARIES "quadmath")
    set(HAVE_QUADMATH "${QUADMATH_FOUND}")
    if (STD_HAVE_ABSQUAD)
      add_definitions(-D STD_HAVE_ABSQUAD)
    endif()
  endif()
endif()

if (USE_QUADMATH AND NOT QUADMATH_FOUND)
  message(FATAL_ERROR "Quadruple precision math support was explicitly requested but is unavailable!")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QUADMATH
  QUADMATH_LIBRARIES HAVE_QUADMATH
)
