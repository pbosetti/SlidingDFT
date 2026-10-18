# R Wrapper

This folder contains an R package that wraps the native `SlidingDFT` library
via [Rcpp](https://www.rcpp.org). The package compiles the C++ sources itself,
so no prior CMake build is needed. A C++20 compiler is required (Xcode command
line tools on macOS, Rtools on Windows, a recent GCC or Clang on Linux).

## Install

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("pbosetti/SlidingDFT", subdir = "r")
```

From a local clone (repository root):

```sh
R CMD INSTALL r
```

## Use

```r
library(SlidingDFT)

window_size <- 256
sdft <- SlidingDFT$new(window_size)        # or SlidingDFT$new(256, 4096)

sdft$update(sin(2 * pi * 17 * (0:4999) / window_size))
power <- sdft$get_power()

sdft$bins                                  # 129
which.max(power) - 1                       # 17
```

`$update()` accepts a numeric vector of any length (a single sample is a
vector of length one); the loop over samples runs in C++.

Notes:

- The R package always uses double precision.
- SIMD acceleration is used on ARM64 (NEON). On x86_64 the package is built
  with R's default compiler flags, which do not enable AVX2, so the scalar
  code path is used.
- Objects wrap native memory and cannot be saved with `saveRDS()` or restored
  across sessions; recreate them instead.

## Development

The library sources under `inst/include/SlidingDFT/` and `src/sliding_dft.cpp`
are copies of the top-level `include/` and `src/` files, because `R CMD build`
only ships files inside the package directory. After changing the core
library, resync from the repository root:

```sh
Rscript r/tools/sync_sources.R
```

CI fails if the copies are out of date. To regenerate the Rcpp glue code and
documentation after changing the wrapper:

```sh
cd r
Rscript -e 'Rcpp::compileAttributes(); roxygen2::roxygenise()'
```

Run the checks:

```sh
R CMD build r
R CMD check SlidingDFT_*.tar.gz
```
