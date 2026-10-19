cmake_minimum_required(VERSION 3.20)

project(SlidingDFT LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(CTest)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(SLIDINGDFT_USE_DOUBLE "Build SlidingDFT with double-precision samples" OFF)
option(SLIDINGDFT_BUILD_BENCHMARKS "Build benchmarks against pocketfft" OFF)

add_library(SlidingDFT STATIC
  src/sliding_dft.cpp
)

add_library(SlidingDFT::SlidingDFT ALIAS SlidingDFT)

add_library(SlidingDFTC SHARED
  src/c_api.cpp
)

add_library(SlidingDFT::SlidingDFTC ALIAS SlidingDFTC)

target_include_directories(SlidingDFT
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(SlidingDFTC
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(SlidingDFT
  PUBLIC
    $<$<BOOL:${SLIDINGDFT_USE_DOUBLE}>:SLIDINGDFT_USE_DOUBLE>
)

target_compile_definitions(SlidingDFTC
  PUBLIC
    $<$<BOOL:${SLIDINGDFT_USE_DOUBLE}>:SLIDINGDFT_USE_DOUBLE>
    SLIDINGDFT_C_EXPORTS
)

target_link_libraries(SlidingDFTC PRIVATE SlidingDFT::SlidingDFT)

add_executable(sliding_dft_example
  src/example.cpp
)

target_link_libraries(sliding_dft_example PRIVATE SlidingDFT::SlidingDFT)

if(SLIDINGDFT_BUILD_BENCHMARKS)
  include(FetchContent)

  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.5
  )

  # pocketfft cpp branch
  FetchContent_Declare(
    pocketfft
    GIT_REPOSITORY https://github.com/mreineck/pocketfft.git
    GIT_TAG 4d865cba24359b05b37da4e75aee24498c3b78e7 # branch: cpp
    SOURCE_SUBDIR does_not_exist
  )

  FetchContent_MakeAvailable(googlebenchmark pocketfft)

  find_package(Threads REQUIRED)

  add_library(pocketfft INTERFACE)
  target_include_directories(pocketfft INTERFACE ${pocketfft_SOURCE_DIR})
  target_link_libraries(pocketfft INTERFACE Threads::Threads)

  add_executable(sliding_dft_benchmark
    bench/benchmark.cpp
  )

  target_link_libraries(sliding_dft_benchmark
    PRIVATE
      SlidingDFT::SlidingDFT
      benchmark::benchmark
      pocketfft
  )
endif()

if(BUILD_TESTING)
  add_test(
    NAME sliding_dft_example_runs
    COMMAND sliding_dft_example
  )
endif()

if(MSVC)
  target_compile_options(SlidingDFT PRIVATE /W4 /permissive-)
  target_compile_options(SlidingDFTC PRIVATE /W4 /permissive-)
  target_compile_options(sliding_dft_example PRIVATE /W4 /permissive-)
  if(TARGET sliding_dft_benchmark)
    target_compile_options(sliding_dft_benchmark PRIVATE /W4 /permissive-)
  endif()
else()
  target_compile_options(SlidingDFT PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(SlidingDFTC PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(sliding_dft_example PRIVATE -Wall -Wextra -Wpedantic)
  if(TARGET sliding_dft_benchmark)
    target_compile_options(sliding_dft_benchmark PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|X86_64|amd64|AMD64|i[3-6]86)$")
  target_compile_options(SlidingDFT PRIVATE -mavx2)
  target_compile_options(SlidingDFTC PRIVATE -mavx2)
  target_compile_options(sliding_dft_example PRIVATE -mavx2)
  if(TARGET sliding_dft_benchmark)
    target_compile_options(sliding_dft_benchmark PRIVATE -mavx2)
  endif()
endif()

set_target_properties(SlidingDFT PROPERTIES OUTPUT_NAME SlidingDFT)
set_target_properties(SlidingDFT PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(SlidingDFTC PROPERTIES OUTPUT_NAME SlidingDFTC)

install(
  TARGETS SlidingDFT SlidingDFTC
  EXPORT SlidingDFTTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
  DIRECTORY include/SlidingDFT
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT SlidingDFTTargets
  FILE SlidingDFTTargets.cmake
  NAMESPACE SlidingDFT::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SlidingDFT
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SlidingDFTConfigVersion.cmake"
  VERSION 1.0.0
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SlidingDFTConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/SlidingDFTConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SlidingDFT
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SlidingDFTConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SlidingDFTConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SlidingDFT
)
