#include <SlidingDFT/sliding_dft.hpp>
#include <benchmark/benchmark.h>
#include <pocketfft_hdronly.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <random>
#include <span>
#include <type_traits>
#include <vector>

namespace {

constexpr std::size_t SignalLength = 65536;
constexpr std::size_t SignalMask = SignalLength - 1;

const std::vector<data_t> &signal() {
  static const std::vector<data_t> samples = [] {
    std::vector<data_t> generated(SignalLength);
    std::mt19937 gen(42);
    std::normal_distribution<data_t> dist(data_t(0), data_t(1));
    for (auto &x : generated) {
      x = dist(gen);
    }
    return generated;
  }();
  return samples;
}

// Sliding-window FFT baseline: keeps the last N samples in a ring buffer and
// computes the full real FFT of the window on demand with pocketfft.
class FftStream {
public:
  explicit FftStream(std::size_t n)
      : _n(n), _ring(n, data_t(0)), _work(n, data_t(0)), _plan(n) {}

  void push(data_t sample) {
    _ring[_index] = sample;
    _index = (_index + 1 == _n) ? 0 : _index + 1;
  }

  // Power spectrum |X_k|^2, k = 0..N/2. Rotation of the window does not
  // change the power spectrum, so the ring buffer is copied as is.
  void get_power(std::vector<data_t> &out) {
    std::copy(_ring.begin(), _ring.end(), _work.begin());
    _plan.exec(_work.data(), data_t(1), true);
    out.resize(_n / 2 + 1);
    // FFTPACK halfcomplex layout: r0, r1, i1, r2, i2, ..., [r_{N/2}]
    out[0] = _work[0] * _work[0];
    for (std::size_t k = 1; k < (_n + 1) / 2; ++k) {
      const data_t re = _work[2 * k - 1];
      const data_t im = _work[2 * k];
      out[k] = re * re + im * im;
    }
    if (_n % 2 == 0) {
      out[_n / 2] = _work[_n - 1] * _work[_n - 1];
    }
  }

private:
  std::size_t _n;
  std::vector<data_t> _ring;
  std::size_t _index = 0;
  std::vector<data_t> _work;
  pocketfft::detail::pocketfft_r<data_t> _plan;
};

void BM_SlidingDFT_Update(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto &sig = signal();
  SlidingDFT sdft(n);
  std::size_t i = 0;
  for (; i < n; ++i) {
    sdft.update(sig[i & SignalMask]);
  }
  for (auto _ : state) {
    sdft.update(sig[i & SignalMask]);
    ++i;
  }
  SlidingDFT::vec_t out;
  sdft.get_power(out);
  benchmark::DoNotOptimize(out.data());
  state.SetItemsProcessed(state.iterations());
}

void BM_SlidingDFT(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto hop = static_cast<std::size_t>(state.range(1));
  const auto &sig = signal();
  SlidingDFT sdft(n);
  SlidingDFT::vec_t out;
  std::size_t i = 0;
  for (; i < n; ++i) {
    sdft.update(sig[i & SignalMask]);
  }
  std::size_t since_power = 0;
  for (auto _ : state) {
    sdft.update(sig[i & SignalMask]);
    ++i;
    if (++since_power == hop) {
      since_power = 0;
      sdft.get_power(out);
      benchmark::DoNotOptimize(out.data());
    }
  }
  state.SetItemsProcessed(state.iterations());
}

void BM_SlidingDFT_Batch(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto hop = static_cast<std::size_t>(state.range(1));
  const auto &sig = signal();
  SlidingDFT sdft(n);
  SlidingDFT::vec_t out;
  sdft.update(std::span(sig.data(), n));
  std::size_t start = n;
  for (auto _ : state) {
    if (start + hop > SignalLength) {
      start = 0;
    }
    sdft.update(std::span(sig.data() + start, hop));
    start += hop;
    sdft.get_power(out);
    benchmark::DoNotOptimize(out.data());
  }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(hop));
}

void BM_PocketFFT(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto hop = static_cast<std::size_t>(state.range(1));
  const auto &sig = signal();
  FftStream fft(n);
  std::vector<data_t> out;
  std::size_t i = 0;
  for (; i < n; ++i) {
    fft.push(sig[i & SignalMask]);
  }
  std::size_t since_power = 0;
  for (auto _ : state) {
    fft.push(sig[i & SignalMask]);
    ++i;
    if (++since_power == hop) {
      since_power = 0;
      fft.get_power(out);
      benchmark::DoNotOptimize(out.data());
    }
  }
  state.SetItemsProcessed(state.iterations());
}

void update_args(benchmark::Benchmark *b) {
  for (const std::int64_t n : {64, 256, 300, 1000, 1024, 4096}) {
    b->Arg(n);
  }
}

void window_and_hop_args(benchmark::Benchmark *b) {
  for (const std::int64_t n : {64, 256, 300, 1000, 1024, 4096}) {
    b->Args({n, 1});
    b->Args({n, 16});
    if (n / 4 != 16) {
      b->Args({n, n / 4});
    }
  }
}

bool check_agreement() {
  constexpr std::size_t Sizes[] = {256, 300};
  constexpr std::size_t Samples = 10000;
  constexpr double Tolerance = std::is_same_v<data_t, float> ? 1e-3 : 1e-9;
  const auto &sig = signal();
  bool ok = true;
  for (const std::size_t n : Sizes) {
    SlidingDFT sdft(n);
    FftStream fft(n);
    for (std::size_t i = 0; i < Samples; ++i) {
      sdft.update(sig[i]);
      fft.push(sig[i]);
    }
    SlidingDFT::vec_t a;
    sdft.get_power(a);
    std::vector<data_t> b;
    fft.get_power(b);

    // Block updates in uneven chunks must match per-sample updates exactly.
    SlidingDFT batch(n);
    constexpr std::size_t Chunks[] = {1, 7, 13};
    for (std::size_t i = 0, c = 0; i < Samples; ++c) {
      const std::size_t len = std::min(Chunks[c % 3], Samples - i);
      batch.update(std::span(sig.data() + i, len));
      i += len;
    }
    SlidingDFT::vec_t a_batch;
    batch.get_power(a_batch);
    const bool identical = a_batch == a;

    double max_diff = 0.0;
    double max_ref = 0.0;
    for (std::size_t k = 0; k < b.size(); ++k) {
      const double ak = static_cast<double>(a[k]);
      const double bk = static_cast<double>(b[k]);
      max_diff = std::max(max_diff, std::fabs(ak - bk));
      max_ref = std::max(max_ref, bk);
    }
    const double rel_err = max_diff / max_ref;
    const bool pass = rel_err <= Tolerance && identical;
    std::fprintf(stderr, "check N=%zu max_rel_err=%.3e tol=%.1e batch=%s %s\n",
                 n, rel_err, Tolerance, identical ? "identical" : "DIFFERENT",
                 pass ? "OK" : "FAIL");
    ok = ok && pass;
  }
  return ok;
}

} // namespace

BENCHMARK(BM_SlidingDFT_Update)->Apply(update_args)->ArgNames({"N"});
BENCHMARK(BM_SlidingDFT)->Apply(window_and_hop_args)->ArgNames({"N", "hop"});
BENCHMARK(BM_SlidingDFT_Batch)
    ->Apply(window_and_hop_args)
    ->ArgNames({"N", "hop"});
BENCHMARK(BM_PocketFFT)->Apply(window_and_hop_args)->ArgNames({"N", "hop"});

int main(int argc, char **argv) {
  if (!check_agreement()) {
    return 1;
  }
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
