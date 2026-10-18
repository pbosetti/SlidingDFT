#' @keywords internal
"_PACKAGE"

#' @useDynLib SlidingDFT, .registration = TRUE
#' @importFrom Rcpp loadModule
#' @importFrom methods new show setMethod evalqOnLoad
NULL

#' Sliding discrete Fourier transform
#'
#' A reference class wrapping the C++ `SlidingDFT` engine. It maintains a
#' running DFT over the most recent `window_size` samples, so the power
#' spectrum can be read after every new sample at a cost proportional to the
#' number of bins.
#'
#' @section Constructor:
#' \describe{
#'   \item{`SlidingDFT$new(window_size, renorm_interval = 4096)`}{
#'     `window_size` is the number of samples in the analysis window (any
#'     positive integer, not necessarily a power of two). `renorm_interval`
#'     is the number of updates between phase renormalizations, which bound
#'     floating-point drift.}
#' }
#'
#' @section Methods and fields:
#' \describe{
#'   \item{`$update(samples)`}{Push a numeric vector of samples, in order.
#'     The loop runs in C++, so passing long vectors is efficient.}
#'   \item{`$get_power()`}{Return the current power spectrum as a numeric
#'     vector of length `bins`. Bin `k` (1-based) is frequency `(k - 1) / N`
#'     cycles per sample.}
#'   \item{`$bins`}{Read-only number of bins, `window_size \%/\% 2 + 1`.}
#' }
#'
#' @details Computations run in double precision. Objects hold a pointer to
#'   native memory and cannot be serialized: after `saveRDS()`/`readRDS()` or
#'   a session restart they must be recreated.
#'
#' @name SlidingDFT
#' @aliases SlidingDFT Rcpp_SlidingDFT-class
#' @export SlidingDFT
#' @examples
#' n <- 256
#' sdft <- SlidingDFT$new(n)
#' sdft$update(sin(2 * pi * 17 * (0:4999) / n))
#' power <- sdft$get_power()
#' which.max(power) - 1 # peak at bin 17
NULL

Rcpp::loadModule("sliding_dft_module", TRUE)

evalqOnLoad({
  setMethod("show", "Rcpp_SlidingDFT", function(object) {
    cat("<SlidingDFT> bins:", object$bins, "\n")
  })
})
