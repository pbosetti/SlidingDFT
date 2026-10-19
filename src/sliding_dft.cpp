#include <SlidingDFT/sliding_dft.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

// MSVC never defines __FMA__, but its /arch:AVX2 implies FMA3 support.
#if defined(__AVX2__) &&                                                       \
    (defined(__FMA__) || (defined(_MSC_VER) && !defined(__clang__)))
#define SLIDINGDFT_HAS_FMA
#endif

namespace {

// Minimal SIMD layer, one backend per target. mul_add(a, b, c) is a * b + c and
// neg_mul_add(a, b, c) is c - a * b.
namespace simd {

#if defined(__AVX2__) && defined(SLIDINGDFT_USE_DOUBLE)
using vec = __m256d;
constexpr std::size_t width = VEC_WIDTH;
inline vec load(const data_t *p) { return _mm256_load_pd(p); }
inline void store(data_t *p, vec v) { _mm256_store_pd(p, v); }
inline vec splat(data_t x) { return _mm256_set1_pd(x); }
inline vec add(vec a, vec b) { return _mm256_add_pd(a, b); }
inline vec mul(vec a, vec b) { return _mm256_mul_pd(a, b); }
#if defined(SLIDINGDFT_HAS_FMA)
inline vec mul_add(vec a, vec b, vec c) { return _mm256_fmadd_pd(a, b, c); }
inline vec neg_mul_add(vec a, vec b, vec c) {
  return _mm256_fnmadd_pd(a, b, c);
}
#else
inline vec mul_add(vec a, vec b, vec c) { return add(mul(a, b), c); }
inline vec neg_mul_add(vec a, vec b, vec c) {
  return _mm256_sub_pd(c, mul(a, b));
}
#endif
#elif defined(__AVX2__)
using vec = __m256;
constexpr std::size_t width = VEC_WIDTH;
inline vec load(const data_t *p) { return _mm256_load_ps(p); }
inline void store(data_t *p, vec v) { _mm256_store_ps(p, v); }
inline vec splat(data_t x) { return _mm256_set1_ps(x); }
inline vec add(vec a, vec b) { return _mm256_add_ps(a, b); }
inline vec mul(vec a, vec b) { return _mm256_mul_ps(a, b); }
#if defined(SLIDINGDFT_HAS_FMA)
inline vec mul_add(vec a, vec b, vec c) { return _mm256_fmadd_ps(a, b, c); }
inline vec neg_mul_add(vec a, vec b, vec c) {
  return _mm256_fnmadd_ps(a, b, c);
}
#else
inline vec mul_add(vec a, vec b, vec c) { return add(mul(a, b), c); }
inline vec neg_mul_add(vec a, vec b, vec c) {
  return _mm256_sub_ps(c, mul(a, b));
}
#endif
#elif defined(__ARM_NEON) && defined(SLIDINGDFT_USE_DOUBLE)
using vec = float64x2_t;
constexpr std::size_t width = VEC_WIDTH;
inline vec load(const data_t *p) { return vld1q_f64(p); }
inline void store(data_t *p, vec v) { vst1q_f64(p, v); }
inline vec splat(data_t x) { return vdupq_n_f64(x); }
inline vec add(vec a, vec b) { return vaddq_f64(a, b); }
inline vec mul(vec a, vec b) { return vmulq_f64(a, b); }
inline vec mul_add(vec a, vec b, vec c) { return vfmaq_f64(c, a, b); }
inline vec neg_mul_add(vec a, vec b, vec c) { return vfmsq_f64(c, a, b); }
#elif defined(__ARM_NEON)
using vec = float32x4_t;
constexpr std::size_t width = VEC_WIDTH;
inline vec load(const data_t *p) { return vld1q_f32(p); }
inline void store(data_t *p, vec v) { vst1q_f32(p, v); }
inline vec splat(data_t x) { return vdupq_n_f32(x); }
inline vec add(vec a, vec b) { return vaddq_f32(a, b); }
inline vec mul(vec a, vec b) { return vmulq_f32(a, b); }
inline vec mul_add(vec a, vec b, vec c) { return vfmaq_f32(c, a, b); }
inline vec neg_mul_add(vec a, vec b, vec c) { return vfmsq_f32(c, a, b); }
#else
using vec = data_t;
constexpr std::size_t width = 1;
inline vec load(const data_t *p) { return *p; }
inline void store(data_t *p, vec v) { *p = v; }
inline vec splat(data_t x) { return x; }
inline vec add(vec a, vec b) { return a + b; }
inline vec mul(vec a, vec b) { return a * b; }
inline vec mul_add(vec a, vec b, vec c) { return a * b + c; }
inline vec neg_mul_add(vec a, vec b, vec c) { return c - a * b; }
#endif

} // namespace simd

using Complex = std::complex<double>;

// Plain complex product: std::complex operator* goes through the C99 Annex G
// NaN/Inf recovery path (__muldc3), which is much slower.
inline Complex mul(Complex a, Complex b) {
  return {a.real() * b.real() - a.imag() * b.imag(),
          a.real() * b.imag() + a.imag() * b.real()};
}

// In-place forward FFT of power-of-two size m, twiddle[j] = e^{-2 pi i j / m}.
void fft_radix2(Complex *x, std::size_t m, const Complex *twiddle) {
  for (std::size_t i = 1, j = 0; i < m; ++i) {
    std::size_t bit = m >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(x[i], x[j]);
    }
  }
  for (std::size_t len = 2; len <= m; len <<= 1) {
    const std::size_t half = len / 2;
    const std::size_t step = m / len;
    for (std::size_t i = 0; i < m; i += len) {
      for (std::size_t j = 0; j < half; ++j) {
        const Complex t = mul(x[i + j + half], twiddle[j * step]);
        x[i + j + half] = x[i + j] - t;
        x[i + j] += t;
      }
    }
  }
}

// State arrays are padded to whole SIMD vectors, so the per-sample kernel has
// no scalar tail. Padded bins have w = 0 and stay exactly zero.
constexpr std::size_t padded_bins(std::size_t bins) {
  return (bins + simd::width - 1) / simd::width * simd::width;
}

// Count consecutive updates y <- (y + delta_j) * w of one bin vector, keeping y
// in registers, with w = c - i s. Every caller goes through this kernel, so
// pushing samples one at a time or in blocks gives identical results.
template <std::size_t Count>
inline void rotate_bins(data_t *yr, data_t *yi, const data_t *wc,
                        const data_t *ws, const simd::vec *deltas) {
  simd::vec r = simd::load(yr);
  simd::vec i = simd::load(yi);
  const simd::vec c = simd::load(wc);
  const simd::vec s = simd::load(ws);
  for (std::size_t j = 0; j < Count; ++j) {
    const simd::vec a = simd::add(r, deltas[j]);
    r = simd::mul_add(a, c, simd::mul(i, s));     // a c + i s
    i = simd::neg_mul_add(a, s, simd::mul(i, c)); // i c - a s
  }
  simd::store(yr, r);
  simd::store(yi, i);
}

} // namespace

SlidingDFT::SlidingDFT(std::size_t window_size, std::size_t renorm_interval)
    : _n(window_size), _k(window_size / 2 + 1), _buffer(window_size, 0),
      _yr(padded_bins(_k), 0), _yi(padded_bins(_k), 0),
      _cos(padded_bins(_k), 0), _sin(padded_bins(_k), 0), _index(0),
      _until_renorm(renorm_interval), _renorm_interval(renorm_interval) {
  if (_n == 0) {
    throw std::invalid_argument("N > 0 required");
  }
  if (_renorm_interval == 0) {
    throw std::invalid_argument("renorm_interval > 0 required");
  }

  // Twiddles w_k = e^{-i 2 pi k / N}, computed in double and then rounded.
  for (std::size_t k = 0; k < _k; ++k) {
    const double angle = 2.0 * std::numbers::pi_v<double> *
                         static_cast<double>(k) / static_cast<double>(_n);
    _cos[k] = static_cast<data_t>(std::cos(angle));
    _sin[k] = static_cast<data_t>(std::sin(angle));
  }

  _fft_size = std::has_single_bit(_n) ? _n : std::bit_ceil(2 * _n - 1);
  _fft_twiddle.resize(_fft_size / 2);
  for (std::size_t j = 0; j < _fft_twiddle.size(); ++j) {
    _fft_twiddle[j] = std::polar(1.0, -2.0 * std::numbers::pi_v<double> *
                                          static_cast<double>(j) /
                                          static_cast<double>(_fft_size));
  }
  _fft_work.resize(_fft_size);

  if (_fft_size != _n) {
    // Bluestein: chirp c_n = e^{i pi n^2 / N}, with n^2 reduced mod 2N to keep
    // the angle accurate, and the FFT of the zero-padded symmetric chirp.
    _chirp.resize(_n);
    for (std::size_t n = 0; n < _n; ++n) {
      const auto n2 = (static_cast<unsigned long long>(n) * n) % (2 * _n);
      _chirp[n] =
          std::polar(1.0, std::numbers::pi_v<double> * static_cast<double>(n2) /
                              static_cast<double>(_n));
    }
    _chirp_fft.assign(_fft_size, 0.0);
    _chirp_fft[0] = _chirp[0];
    for (std::size_t n = 1; n < _n; ++n) {
      _chirp_fft[n] = _chirp[n];
      _chirp_fft[_fft_size - n] = _chirp[n];
    }
    fft_radix2(_chirp_fft.data(), _fft_size, _fft_twiddle.data());
  }
}

void SlidingDFT::update(data_t sample) {
  const data_t delta = sample - _buffer[_index];
  _buffer[_index] = sample;
  if (++_index == _n) {
    _index = 0;
  }

  update_bins<1>(&delta);

  if (--_until_renorm == 0) {
    _until_renorm = _renorm_interval;
    recompute_fft_like();
  }
}

void SlidingDFT::update(std::span<const data_t> samples) {
  // Blocks of samples share one pass over the bins, with the bins held in
  // registers. A block never straddles a recomputation, so the result is the
  // same as with per-sample updates. Larger blocks run out of registers.
  constexpr std::size_t Block = 4;
  data_t deltas[Block] = {};
  std::size_t pos = 0;
  while (pos < samples.size()) {
    const std::size_t count =
        std::min({Block, samples.size() - pos, _until_renorm});
    for (std::size_t j = 0; j < count; ++j) {
      const data_t sample = samples[pos + j];
      deltas[j] = sample - _buffer[_index];
      _buffer[_index] = sample;
      if (++_index == _n) {
        _index = 0;
      }
    }

    if (count == Block) {
      update_bins<Block>(deltas);
    } else {
      for (std::size_t j = 0; j < count; ++j) {
        update_bins<1>(deltas + j);
      }
    }

    pos += count;
    _until_renorm -= count;
    if (_until_renorm == 0) {
      _until_renorm = _renorm_interval;
      recompute_fft_like();
    }
  }
}

void SlidingDFT::get_power(vec_t &power_spectrum) const {
  power_spectrum.resize(_k);
  power(power_spectrum.data());
}

std::size_t SlidingDFT::bins() const noexcept { return _k; }

// Kernels work on local copies of the member pointers and sizes: SIMD stores
// may alias any object, so going through the members would reload them from
// memory on every iteration.

template <std::size_t Count>
void SlidingDFT::update_bins(const data_t *deltas) {
  data_t *yr = _yr.data();
  data_t *yi = _yi.data();
  const data_t *wc = _cos.data();
  const data_t *ws = _sin.data();
  const std::size_t bins = _yr.size(); // a multiple of simd::width

  simd::vec vdeltas[Count];
  for (std::size_t j = 0; j < Count; ++j) {
    vdeltas[j] = simd::splat(deltas[j]);
  }

  // Two vectors per iteration, then at most one more.
  constexpr std::size_t w = simd::width;
  std::size_t k = 0;
  for (; k + 2 * w <= bins; k += 2 * w) {
    rotate_bins<Count>(yr + k, yi + k, wc + k, ws + k, vdeltas);
    rotate_bins<Count>(yr + k + w, yi + k + w, wc + k + w, ws + k + w, vdeltas);
  }
  if (k < bins) {
    rotate_bins<Count>(yr + k, yi + k, wc + k, ws + k, vdeltas);
  }
}

void SlidingDFT::power(data_t *out) const {
  const data_t *yr = _yr.data();
  const data_t *yi = _yi.data();
  const std::size_t bins = _k;

  // |X_k|^2 = |y_k|^2: y_k is the conjugate of the DFT of the window ordered
  // from the oldest sample, so no phase reference is needed.
  const auto norm = [&](std::size_t k) {
    const simd::vec r = simd::load(yr + k);
    const simd::vec i = simd::load(yi + k);
    simd::store(out + k, simd::mul_add(r, r, simd::mul(i, i)));
  };
  std::size_t k = 0;
  for (; k + 2 * simd::width <= bins; k += 2 * simd::width) {
    norm(k);
    norm(k + simd::width);
  }
  if (k + simd::width <= bins) {
    norm(k);
    k += simd::width;
  }
  for (; k < bins; ++k) {
    out[k] = yr[k] * yr[k] + yi[k] * yi[k];
  }
}

void SlidingDFT::recompute_fft_like() {
  // DFT X_k of the window ordered from the oldest sample, in double precision.
  Complex *x = _fft_work.data();
  for (std::size_t n = 0; n < _n; ++n) {
    const double sample = _buffer[(_index + n) % _n];
    x[n] = _chirp.empty() ? Complex(sample) : sample * std::conj(_chirp[n]);
  }

  if (_chirp.empty()) {
    fft_radix2(x, _fft_size, _fft_twiddle.data());
  } else {
    // X_k = conj(c_k) * (a * c)_k, the circular convolution computed as
    // ifft(fft(a) . fft(c)), with ifft(z) = conj(fft(conj(z))) / M.
    std::fill(x + _n, x + _fft_size, Complex(0.0));
    fft_radix2(x, _fft_size, _fft_twiddle.data());
    for (std::size_t m = 0; m < _fft_size; ++m) {
      x[m] = std::conj(mul(x[m], _chirp_fft[m]));
    }
    fft_radix2(x, _fft_size, _fft_twiddle.data());
    const double scale = 1.0 / static_cast<double>(_fft_size);
    for (std::size_t k = 0; k < _k; ++k) {
      x[k] = mul(std::conj(x[k]) * scale, std::conj(_chirp[k]));
    }
  }

  // The recursion in update_bins() tracks the conjugate of X_k.
  for (std::size_t k = 0; k < _k; ++k) {
    _yr[k] = static_cast<data_t>(x[k].real());
    _yi[k] = static_cast<data_t>(-x[k].imag());
  }
}
