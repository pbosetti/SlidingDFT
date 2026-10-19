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
- Objects wrap native memory and cannot be saved with `saveRDS()` or restored
  across sessions; recreate them instead.

## Build Configuration

The package is compiled with R's own C++20 compiler and flags, which you can
inspect with `R CMD config CXX20` and `R CMD config CXX20FLAGS`. Extra flags
come from a user Makevars file, and `+=` appends them to R's defaults.

### ARM64

Nothing to configure: NEON is always available, and `-O3` gives no measurable
gain over R's default `-O2`.

### x86_64 on Linux and macOS: enable AVX2

R's default flags target a generic x86_64 CPU, so the scalar code path is
used. On a CPU with AVX2 and FMA (Intel since Haswell, 2013; AMD since Zen,
2017) the SIMD path is much faster. Choose one of the following:

- One install, without touching any configuration file:

  ```r
  withr::with_makevars(
    c(CXX20FLAGS = "-mavx2 -mfma"), assignment = "+=",
    remotes::install_github("pbosetti/SlidingDFT", subdir = "r", force = TRUE)
  )
  ```

  `force = TRUE` rebuilds the package even if the same version is already
  installed.

- Permanently, for every package you build, in `~/.R/Makevars`:

  ```make
  CXX20FLAGS += -mavx2 -mfma
  ```

- From the shell, with a Makevars file used only for this build:

  ```sh
  echo 'CXX20FLAGS += -mavx2 -mfma' > sdft.mk
  R_MAKEVARS_USER="$PWD/sdft.mk" R CMD INSTALL r
  ```

The flags appear on the compiler command lines in the install log
(`install.packages(..., quiet = FALSE)`, or the output of `R CMD INSTALL`).
`-march=native` works too, but ties the build to the build machine's CPU. The
resulting library stops with an illegal-instruction error on a CPU without
AVX2, so don't use these flags for binaries you share.

### Windows (Rtools)

Keep R's default flags. With AVX enabled, the GCC compiler in Rtools can place
32-byte vectors on a stack that is only 16-byte aligned
([GCC bug 54412](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54412)), which
crashes at run time. If you need other flags, R on Windows reads the first file
found among `~/.R/Makevars.ucrt`, `~/.R/Makevars.win64`, `~/.R/Makevars.win`
and `~/.R/Makevars`, unless `R_MAKEVARS_USER` points elsewhere.

### Compiler and precision

To use another compiler, set it in the user Makevars file, for example
`CXX20 = clang++`. Precision cannot be configured: R numbers are doubles, so
`src/Makevars` always defines `SLIDINGDFT_USE_DOUBLE`, and the build fails if
it is removed.

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
