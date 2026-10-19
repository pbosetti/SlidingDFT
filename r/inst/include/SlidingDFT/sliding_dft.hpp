#pragma once

#include <SlidingDFT/allocator.hpp>
#include <SlidingDFT/types.hpp>

#include <complex>
#include <span>

/**
 * @class SlidingDFT
 * @brief Sliding discrete Fourier transform for real-valued input streams.
 *
 * The class maintains a running DFT over a fixed-size window of the most recent
 * samples. SIMD acceleration is selected at compile time when AVX2 or NEON is
 * available; otherwise a scalar implementation is used.
 */
class SlidingDFT {
public:
  using vec_t = std::vector<data_t, AlignedAllocator<data_t, 32>>;

  /**
   * @brief Construct a sliding DFT for a fixed window length.
   * @param window_size Number of samples in the analysis window.
   * @param renorm_interval Number of updates between exact recomputations
   *        of the bins, which bound floating-point drift.
   * @throws std::invalid_argument If @p window_size or @p renorm_interval is
   *         zero.
   */
  SlidingDFT(std::size_t window_size, std::size_t renorm_interval = 4096);

  /**
   * @brief Push one new sample into the window and update all bins.
   * @param sample New sample value.
   */
  void update(data_t sample);

  /**
   * @brief Push a block of samples into the window, oldest first.
   *
   * Gives exactly the same result as calling update() on each sample, but is
   * faster: consecutive samples are applied to each bin in one pass.
   * @param samples New sample values.
   */
  void update(std::span<const data_t> samples);

  /**
   * @brief Compute the current power spectrum.
   * @param power_spectrum Output vector resized to `bins()` elements.
   */
  void get_power(vec_t &power_spectrum) const;

  /**
   * @brief Return the number of frequency bins.
   * @return `window_size / 2 + 1`.
   */
  [[nodiscard]] std::size_t bins() const noexcept;

private:
  template <std::size_t Count> void update_bins(const data_t *deltas);
  void power(data_t *out) const;
  void recompute_fft_like();

  std::size_t _n;
  std::size_t _k;
  vec_t _buffer;
  vec_t _yr;
  vec_t _yi;
  // Twiddles w_k = e^{-i 2 pi k / N} = _cos[k] - i _sin[k].
  vec_t _cos;
  vec_t _sin;
  std::size_t _index;
  std::size_t _until_renorm;
  std::size_t _renorm_interval;
  // FFT used by recompute_fft_like(): radix-2 of size _fft_size, with
  // Bluestein's algorithm when the window size is not a power of two.
  std::size_t _fft_size;
  std::vector<std::complex<double>> _fft_twiddle;
  std::vector<std::complex<double>> _chirp;
  std::vector<std::complex<double>> _chirp_fft;
  std::vector<std::complex<double>> _fft_work;
};

using SlidingDFT_AVX2 = SlidingDFT;
