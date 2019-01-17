add_executable(fdtd main.c fdtd.c fdtd1D.c fdtd2D.c fdtd3D.c initialize.c fdtd_common.c)
target_include_directories(fdtd PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_definitions(fdtd PRIVATE -DFDTD_USE_DOUBLE)
target_link_libraries(fdtd PRIVATE m)
set_property(TARGET fdtd
             PROPERTY C_STANDARD 11)

find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(fdtd PRIVATE OpenMP::OpenMP_C)
endif()

# Compile Options
include(compile-flags-helpers)
include(${PROJECT_SOURCE_DIR}/optimization_flags.cmake)

if (DEFINED ADDITIONAL_BENCHMARK_COMPILE_OPTIONS)
  add_compiler_option_to_target_type(fdtd Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_COMPILE_OPTIONS})
endif()

foreach(compile_type IN ITEMS Release RelWithDebInfo)
  add_compiler_option_to_target_type(fdtd ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_COMPILE_OPTIONS})
  add_linker_option_to_target_type(fdtd ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_LINK_OPTIONS})
endforeach()

add_compiler_option_to_target_type(fdtd Debug PRIVATE ${ADDITIONAL_DEBUG_COMPILE_OPTIONS})

# Linker Options

if (DEFINED ADDITIONAL_BENCHMARK_LINK_OPTIONS)
  add_linker_option_to_target_type(fdtd Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_LINK_OPTIONS})
endif()

add_sanitizers_to_target(fdtd Debug PRIVATE address undefined)

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if((result) AND USE_IPO)
  set_property(TARGET fdtd PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

install(TARGETS fdtd RUNTIME DESTINATION bin)
