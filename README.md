# SlidingDFT

[![CI](https://github.com/pbosetti/SlidingDFT/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/pbosetti/SlidingDFT/actions/workflows/main.yml)
[![R](https://github.com/pbosetti/SlidingDFT/actions/workflows/r.yml/badge.svg?branch=main)](https://github.com/pbosetti/SlidingDFT/actions/workflows/r.yml)

`SlidingDFT` is a small C++20 library that maintains a sliding discrete Fourier
transform over a fixed-size window of real-valued samples.

The project builds a static library and installs public headers under
`SlidingDFT/`, so consumers can include:

```cpp
#include <SlidingDFT/sliding_dft.hpp>
```

## Features

- C++20 library with CMake build system
- C interface
- Static library target: `SlidingDFT::SlidingDFT`
- SIMD backends where available
- AVX2 on x86_64
- NEON on ARM64, including Apple Silicon
- `float` and `double` sample precision via CMake option
- Installable CMake package and `FetchContent`-friendly layout
- Python interface (see [python/README.md](python/README.md))
- R interface via Rcpp (see [r/README.md](r/README.md))
- Benchmark against pocketfft (optional)

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

## Benchmark

The benchmark is optional and off by default. It fetches Google Benchmark and
[pocketfft](https://github.com/mreineck/pocketfft) with `FetchContent`:

```sh
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DSLIDINGDFT_BUILD_BENCHMARKS=ON
cmake --build build
```

Run all benchmarks, or filter by window size or benchmark name:

```sh
./build/sliding_dft_benchmark
./build/sliding_dft_benchmark --benchmark_filter='N:256/'
./build/sliding_dft_benchmark --benchmark_filter='BM_PocketFFT'
```

There are three benchmarks, all reporting input samples per second:

- `BM_SlidingDFT_Update/N`: only `update()` for each sample.
- `BM_SlidingDFT/N/hop`: `update()` for each sample, and `get_power()` every
  `hop` samples.
- `BM_PocketFFT/N/hop`: a ring buffer write for each sample and, every `hop`
  samples, a real FFT of the last `N` samples followed by the power
  computation. It produces the same output as `BM_SlidingDFT`.

`N` ranges over 64 to 4096, including non-powers of two, and `hop` is 1, 16,
or `N/4`.

Before timing, the executable checks that both implementations produce the
same spectrum, and exits with a non-zero status if they do not.

To compare the two, look at the samples per second of `BM_SlidingDFT` and
`BM_PocketFFT` at the same `N` and `hop`. A sliding DFT costs O(N) per sample,
while an FFT costs O(N log N) per transform, so the sliding DFT should win at
small hops. At large hops the cost of the FFT is amortized over many samples,
so the FFT can win.

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
