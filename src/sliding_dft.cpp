#include <SlidingDFT/sliding_dft.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace {

#if defined(__AVX2__) && defined(SLIDINGDFT_USE_DOUBLE)
using vec = __m256d;

inline vec load(const data_t *p) { return _mm256_load_pd(p); }

inline void store(data_t *p, vec v) { _mm256_store_pd(p, v); }
#elif defined(__AVX2__)
using vec = __m256;

inline vec load(const data_t *p) { return _mm256_load_ps(p); }

inline void store(data_t *p, vec v) { _mm256_store_ps(p, v); }
#elif defined(__ARM_NEON) && defined(SLIDINGDFT_USE_DOUBLE)
using vec = float64x2_t;

inline vec load(const data_t *p) { return vld1q_f64(p); }

inline void store(data_t *p, vec v) { vst1q_f64(p, v); }
#elif defined(__ARM_NEON)
using vec = float32x4_t;

inline vec load(const data_t *p) { return vld1q_f32(p); }

inline void store(data_t *p, vec v) { vst1q_f32(p, v); }
#endif

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

} // namespace

SlidingDFT::SlidingDFT(std::size_t window_size, std::size_t renorm_interval)
    : _n(window_size), _k(window_size / 2 + 1), _buffer(window_size, 0),
      _yr(_k, 0), _yi(_k, 0), _wr(_k, 0), _wi(_k, 0), _index(0), _filled(0),
      _count(0), _renorm_interval(renorm_interval) {
  if (_n == 0) {
    throw std::invalid_argument("N > 0 required");
  }

  const data_t two_pi_over_n = static_cast<data_t>(
      2.0 * std::numbers::pi_v<double> / static_cast<double>(_n));

  for (std::size_t k = 0; k < _k; ++k) {
    const data_t angle = two_pi_over_n * static_cast<data_t>(k);
    _wr[k] = std::cos(angle);
    _wi[k] = -std::sin(angle);
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
  const data_t old_sample = _buffer[_index];
  _buffer[_index] = sample;
  _index = (_index + 1) % _n;

  if (_filled < _n) {
    ++_filled;
  }

  const data_t delta = sample - old_sample;
  update_bins(delta);

  if (++_count % _renorm_interval == 0) {
    recompute_fft_like();
  }
}

void SlidingDFT::get_power(vec_t &power_spectrum) const {
  power_spectrum.resize(_k);
  power(power_spectrum.data());
}

std::size_t SlidingDFT::bins() const noexcept { return _k; }

void SlidingDFT::update_bins(data_t delta) {
#if defined(__AVX2__) && defined(SLIDINGDFT_USE_DOUBLE)
  const vec vdelta = _mm256_set1_pd(delta);
  std::size_t k = 0;
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
    const vec wr = load(&_wr[k]);
    const vec wi = load(&_wi[k]);

    yr = _mm256_add_pd(yr, vdelta);

#if defined(__FMA__)
    const vec real = _mm256_fmsub_pd(yr, wr, _mm256_mul_pd(yi, wi));
    const vec imag = _mm256_fmadd_pd(yr, wi, _mm256_mul_pd(yi, wr));
#else
    const vec real =
        _mm256_sub_pd(_mm256_mul_pd(yr, wr), _mm256_mul_pd(yi, wi));
    const vec imag =
        _mm256_add_pd(_mm256_mul_pd(yr, wi), _mm256_mul_pd(yi, wr));
#endif

    store(&_yr[k], real);
    store(&_yi[k], imag);
  }

  for (; k < _k; ++k) {
    const data_t yr = _yr[k] + delta;
    const data_t yi = _yi[k];
    _yr[k] = yr * _wr[k] - yi * _wi[k];
    _yi[k] = yr * _wi[k] + yi * _wr[k];
  }
#elif defined(__AVX2__)
  const vec vdelta = _mm256_set1_ps(delta);
  std::size_t k = 0;
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
    const vec wr = load(&_wr[k]);
    const vec wi = load(&_wi[k]);

    yr = _mm256_add_ps(yr, vdelta);

#if defined(__FMA__)
    const vec real = _mm256_fmsub_ps(yr, wr, _mm256_mul_ps(yi, wi));
    const vec imag = _mm256_fmadd_ps(yr, wi, _mm256_mul_ps(yi, wr));
#else
    const vec real =
        _mm256_sub_ps(_mm256_mul_ps(yr, wr), _mm256_mul_ps(yi, wi));
    const vec imag =
        _mm256_add_ps(_mm256_mul_ps(yr, wi), _mm256_mul_ps(yi, wr));
#endif

    store(&_yr[k], real);
    store(&_yi[k], imag);
  }

  for (; k < _k; ++k) {
    const data_t yr = _yr[k] + delta;
    const data_t yi = _yi[k];
    _yr[k] = yr * _wr[k] - yi * _wi[k];
    _yi[k] = yr * _wi[k] + yi * _wr[k];
  }
#elif defined(__ARM_NEON) && defined(SLIDINGDFT_USE_DOUBLE)
  const vec vdelta = vdupq_n_f64(delta);
  std::size_t k = 0;
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
    const vec wr = load(&_wr[k]);
    const vec wi = load(&_wi[k]);

    yr = vaddq_f64(yr, vdelta);

    const vec real = vsubq_f64(vmulq_f64(yr, wr), vmulq_f64(yi, wi));
    const vec imag = vaddq_f64(vmulq_f64(yr, wi), vmulq_f64(yi, wr));

    store(&_yr[k], real);
    store(&_yi[k], imag);
  }

  for (; k < _k; ++k) {
    const data_t yr = _yr[k] + delta;
    const data_t yi = _yi[k];
    _yr[k] = yr * _wr[k] - yi * _wi[k];
    _yi[k] = yr * _wi[k] + yi * _wr[k];
  }
#elif defined(__ARM_NEON)
  const vec vdelta = vdupq_n_f32(delta);
  std::size_t k = 0;
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
    const vec wr = load(&_wr[k]);
    const vec wi = load(&_wi[k]);

    yr = vaddq_f32(yr, vdelta);

    const vec real = vsubq_f32(vmulq_f32(yr, wr), vmulq_f32(yi, wi));
    const vec imag = vaddq_f32(vmulq_f32(yr, wi), vmulq_f32(yi, wr));

    store(&_yr[k], real);
    store(&_yi[k], imag);
  }

  for (; k < _k; ++k) {
    const data_t yr = _yr[k] + delta;
    const data_t yi = _yi[k];
    _yr[k] = yr * _wr[k] - yi * _wi[k];
    _yi[k] = yr * _wi[k] + yi * _wr[k];
  }
#else
  for (std::size_t k = 0; k < _k; ++k) {
    const data_t yr = _yr[k] + delta;
    const data_t yi = _yi[k];
    _yr[k] = yr * _wr[k] - yi * _wi[k];
    _yi[k] = yr * _wi[k] + yi * _wr[k];
  }
#endif
}

void SlidingDFT::power(data_t *out) const {
  std::size_t k = 0;
#if defined(__AVX2__) && defined(SLIDINGDFT_USE_DOUBLE)
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    const vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
#if defined(__FMA__)
    store(&out[k], _mm256_fmadd_pd(yr, yr, _mm256_mul_pd(yi, yi)));
#else
    store(&out[k], _mm256_add_pd(_mm256_mul_pd(yr, yr), _mm256_mul_pd(yi, yi)));
#endif
  }
#elif defined(__AVX2__)
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    const vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
#if defined(__FMA__)
    store(&out[k], _mm256_fmadd_ps(yr, yr, _mm256_mul_ps(yi, yi)));
#else
    store(&out[k], _mm256_add_ps(_mm256_mul_ps(yr, yr), _mm256_mul_ps(yi, yi)));
#endif
  }
#elif defined(__ARM_NEON) && defined(SLIDINGDFT_USE_DOUBLE)
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    const vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
    store(&out[k], vaddq_f64(vmulq_f64(yr, yr), vmulq_f64(yi, yi)));
  }
#elif defined(__ARM_NEON)
  for (; k + VEC_WIDTH <= _k; k += VEC_WIDTH) {
    const vec yr = load(&_yr[k]);
    const vec yi = load(&_yi[k]);
    store(&out[k], vaddq_f32(vmulq_f32(yr, yr), vmulq_f32(yi, yi)));
  }
#endif
  // |X_k|^2 = |y_k|^2: y_k is the conjugate of the DFT of the window ordered
  // from the oldest sample, so no phase reference is needed.
  for (; k < _k; ++k) {
    out[k] = _yr[k] * _yr[k] + _yi[k] * _yi[k];
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
