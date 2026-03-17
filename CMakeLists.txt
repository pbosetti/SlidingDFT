cmake_minimum_required(VERSION 3.20)

project(SlidingDFT LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SLIDINGDFT_USE_DOUBLE "Build SlidingDFT with double-precision samples" OFF)

add_library(SlidingDFT STATIC
  src/sliding_dft.cpp
)

add_library(SlidingDFT::SlidingDFT ALIAS SlidingDFT)

target_include_directories(SlidingDFT
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(SlidingDFT
  PUBLIC
    $<$<BOOL:${SLIDINGDFT_USE_DOUBLE}>:SLIDINGDFT_USE_DOUBLE>
)

add_executable(sliding_dft_example
  src/example.cpp
)

target_link_libraries(sliding_dft_example PRIVATE SlidingDFT::SlidingDFT)

if(MSVC)
  target_compile_options(SlidingDFT PRIVATE /W4 /permissive-)
  target_compile_options(sliding_dft_example PRIVATE /W4 /permissive-)
else()
  target_compile_options(SlidingDFT PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(sliding_dft_example PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|X86_64|amd64|AMD64|i[3-6]86)$")
  target_compile_options(SlidingDFT PRIVATE -mavx2)
  target_compile_options(sliding_dft_example PRIVATE -mavx2)
endif()

set_target_properties(SlidingDFT PROPERTIES OUTPUT_NAME SlidingDFT)

install(
  TARGETS SlidingDFT
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
