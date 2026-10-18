test_that("bins follow window_size %/% 2 + 1", {
  expect_equal(SlidingDFT$new(256)$bins, 129L)
  expect_equal(SlidingDFT$new(300)$bins, 151L)
  expect_equal(SlidingDFT$new(7, 100)$bins, 4L)
})

test_that("power peaks at the bin of a pure tone", {
  n <- 300
  sdft <- SlidingDFT$new(n)
  sdft$update(sin(2 * pi * 17 * (0:4999) / n))
  power <- sdft$get_power()
  expect_length(power, sdft$bins)
  expect_equal(which.max(power) - 1, 17)
})

test_that("power matches a direct DFT of the window across renormalizations", {
  x <- rnorm(1000)
  for (n in c(64, 101)) {
    for (renorm in c(37, 128, 100000)) {
      sdft <- SlidingDFT$new(n, renorm)
      sdft$update(x)
      ref <- Mod(fft(utils::tail(x, n)))^2
      expect_equal(sdft$get_power(), ref[seq_len(n %/% 2 + 1)], tolerance = 1e-8)
    }
  }
})

test_that("vectorized and per-sample updates agree", {
  x <- rnorm(500)
  a <- SlidingDFT$new(32)
  b <- SlidingDFT$new(32)
  a$update(x)
  for (s in x) b$update(s)
  expect_identical(a$get_power(), b$get_power())
})

test_that("invalid arguments are rejected", {
  expect_error(SlidingDFT$new(0), "window_size")
  expect_error(SlidingDFT$new(-5), "window_size")
  expect_error(SlidingDFT$new(16, 0), "renorm_interval")
})
