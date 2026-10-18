# Copy the core C++ library sources into the R package.
#
# `R CMD build` only ships files that live inside the package directory, so the
# library sources are vendored here. Run from the repository root after any
# change to `include/` or `src/sliding_dft.cpp`:
#
#   Rscript r/tools/sync_sources.R

root <- normalizePath(".")
pkg <- file.path(root, "r")
if (!file.exists(file.path(pkg, "DESCRIPTION"))) {
  stop("run this script from the repository root")
}

headers <- c("allocator.hpp", "sliding_dft.hpp", "types.hpp")
ok <- file.copy(
  file.path(root, "include", "SlidingDFT", headers),
  file.path(pkg, "inst", "include", "SlidingDFT"),
  overwrite = TRUE
)
ok <- c(ok, file.copy(
  file.path(root, "src", "sliding_dft.cpp"),
  file.path(pkg, "src"),
  overwrite = TRUE
))
if (!all(ok)) {
  stop("failed to copy one or more sources")
}
message("SlidingDFT sources synced into r/")
