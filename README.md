# SlidingDFT

`SlidingDFT` is a small C++20 library that maintains a sliding discrete Fourier
transform over a fixed-size window of real-valued samples.

The project builds a static library and installs public headers under
`SlidingDFT/`, so consumers can include:

```cpp
#include <SlidingDFT/sliding_dft.hpp>
```

## Features

- C++20 library with CMake build system
- Static library target: `SlidingDFT::SlidingDFT`
- SIMD backends where available
- AVX2 on x86_64
- NEON on ARM64, including Apple Silicon
- `float` and `double` sample precision via CMake option
- Installable CMake package and `FetchContent`-friendly layout

## Build

Configure and build:

```sh
cmake -B build
cmake --build build
```

Run the example:

```sh
./build/sliding_dft_example
```

Run the smoke test:

```sh
cmake -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Precision

By default, the library uses:

```cpp
using data_t = float;
```

Enable double precision with:

```sh
cmake -B build -DSLIDINGDFT_USE_DOUBLE=ON
cmake --build build
```

This changes `data_t` to `double`. The SIMD implementation supports both
precisions.

## Install

Install to a prefix:

```sh
cmake -B build
cmake --build build
cmake --install build --prefix /path/to/prefix
```

Installed artifacts:

- `lib/libSlidingDFT.a`
- `include/SlidingDFT/*.hpp`
- `lib/cmake/SlidingDFT/*`

## Use With CMake

After installation:

```cmake
find_package(SlidingDFT CONFIG REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE SlidingDFT::SlidingDFT)
```

Example source:

```cpp
#include <SlidingDFT/sliding_dft.hpp>

int main() {
  SlidingDFT sdft(256);
  sdft.update(1.0f);

  SlidingDFT::vec_t power;
  sdft.get_power(power);
  return 0;
}
```

## Use With FetchContent

```cmake
include(FetchContent)

FetchContent_Declare(
  SlidingDFT
  GIT_REPOSITORY https://github.com/your-org/SlidingDFT.git
  GIT_TAG main
)

FetchContent_MakeAvailable(SlidingDFT)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE SlidingDFT::SlidingDFT)
```

If you want the dependency built with double precision, set the option before
`FetchContent_MakeAvailable`:

```cmake
set(SLIDINGDFT_USE_DOUBLE ON CACHE BOOL "" FORCE)
```

## CI

GitHub Actions builds and runs the smoke test on pushes to `main` for:

- Linux x86_64
- Linux arm64
- macOS x86_64
- macOS arm64
