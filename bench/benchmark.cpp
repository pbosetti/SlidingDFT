#include <SlidingDFT/sliding_dft.hpp>
#include <benchmark/benchmark.h>
#include <pocketfft_hdronly.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <random>
#include <span>
#include <string>
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
  const auto power_every = static_cast<std::size_t>(state.range(1));
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
    if (++since_power == power_every) {
      since_power = 0;
      sdft.get_power(out);
      benchmark::DoNotOptimize(out.data());
    }
  }
  state.SetItemsProcessed(state.iterations());
}

// A spectrum every power_every samples, which arrive in update() calls of
// batch samples (the last call of each period may be shorter).
void BM_SlidingDFT_Batch(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto power_every = static_cast<std::size_t>(state.range(1));
  const auto batch = static_cast<std::size_t>(state.range(2));
  const auto &sig = signal();
  SlidingDFT sdft(n);
  SlidingDFT::vec_t out;
  sdft.update(std::span(sig.data(), n));
  std::size_t start = n;
  for (auto _ : state) {
    if (start + power_every > SignalLength) {
      start = 0;
    }
    for (std::size_t pushed = 0; pushed < power_every; pushed += batch) {
      const std::size_t len = std::min(batch, power_every - pushed);
      sdft.update(std::span(sig.data() + start + pushed, len));
    }
    start += power_every;
    sdft.get_power(out);
    benchmark::DoNotOptimize(out.data());
  }
  state.SetItemsProcessed(state.iterations() *
                          static_cast<std::int64_t>(power_every));
}

void BM_PocketFFT(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const auto power_every = static_cast<std::size_t>(state.range(1));
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
    if (++since_power == power_every) {
      since_power = 0;
      fft.get_power(out);
      benchmark::DoNotOptimize(out.data());
    }
  }
  state.SetItemsProcessed(state.iterations());
}

constexpr std::int64_t WindowSizes[] = {64, 256, 300, 1000, 1024, 4096};

void update_args(benchmark::Benchmark *b) {
  for (const std::int64_t n : WindowSizes) {
    b->Arg(n);
  }
}

// Spectrum intervals: every sample, every 16, every N/4 and every N samples
// (non-overlapping windows), without duplicates.
std::vector<std::int64_t> power_intervals(std::int64_t n) {
  std::vector<std::int64_t> out;
  for (const std::int64_t every :
       {std::int64_t{1}, std::int64_t{16}, n / 4, n}) {
    if (every >= 1 && std::find(out.begin(), out.end(), every) == out.end()) {
      out.push_back(every);
    }
  }
  return out;
}

void window_and_power_args(benchmark::Benchmark *b) {
  for (const std::int64_t n : WindowSizes) {
    for (const std::int64_t every : power_intervals(n)) {
      b->Args({n, every});
    }
  }
}

// Batch sizes up to the spectrum interval, which is the largest useful one.
void window_power_and_batch_args(benchmark::Benchmark *b) {
  for (const std::int64_t n : WindowSizes) {
    for (const std::int64_t every : power_intervals(n)) {
      for (const std::int64_t batch : {std::int64_t{1}, std::int64_t{4},
                                       std::int64_t{16}, std::int64_t{64}}) {
        if (batch < every) {
          b->Args({n, every, batch});
        }
      }
      b->Args({n, every, every});
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

// Tidy CSV output for --csv: one row per run, with the benchmark name split
// into variables. batch is the number of samples pushed per update call;
// power_every is the number of samples between spectra (NA: never computed).
class TidyCsvReporter : public benchmark::BenchmarkReporter {
public:
  bool ReportContext(const Context &) override {
    GetOutputStream() << "algorithm,precision,N,batch,power_every,repetition,"
                         "iterations,samples_per_s,ns_per_sample\n";
    return true;
  }

  void ReportRuns(const std::vector<Run> &runs) override {
    for (const Run &run : runs) {
      if (run.run_type != Run::RT_Iteration || run.skipped) {
        continue;
      }
      const auto args = parse_args(run.run_name.args);
      const std::string &fn = run.run_name.function_name;
      const std::string n = args.at("N");
      std::string algorithm = "pocketfft";
      std::string batch = "1";
      std::string power_every = "NA";
      if (fn == "BM_SlidingDFT_Update") {
        algorithm = "sliding_dft";
      } else if (fn == "BM_SlidingDFT") {
        algorithm = "sliding_dft";
        power_every = args.at("power_every");
      } else if (fn == "BM_SlidingDFT_Batch") {
        algorithm = "sliding_dft_batch";
        batch = args.at("batch");
        power_every = args.at("power_every");
      } else {
        power_every = args.at("power_every");
      }
      const double rate = run.counters.at("items_per_second");
      GetOutputStream() << algorithm << ','
                        << (std::is_same_v<data_t, float> ? "float" : "double")
                        << ',' << n << ',' << batch << ',' << power_every << ','
                        << std::max<std::int64_t>(run.repetition_index, 0)
                        << ',' << run.iterations << ',' << rate << ','
                        << 1e9 / rate << '\n';
    }
  }

private:
  // "N:256/power_every:16" -> {N: 256, power_every: 16}
  static std::map<std::string, std::string> parse_args(const std::string &s) {
    std::map<std::string, std::string> out;
    std::size_t start = 0;
    while (start < s.size()) {
      std::size_t end = s.find('/', start);
      if (end == std::string::npos) {
        end = s.size();
      }
      const std::string field = s.substr(start, end - start);
      const std::size_t colon = field.find(':');
      if (colon != std::string::npos) {
        out[field.substr(0, colon)] = field.substr(colon + 1);
      }
      start = end + 1;
    }
    return out;
  }
};

} // namespace

BENCHMARK(BM_SlidingDFT_Update)->Apply(update_args)->ArgNames({"N"});
BENCHMARK(BM_SlidingDFT)
    ->Apply(window_and_power_args)
    ->ArgNames({"N", "power_every"});
BENCHMARK(BM_SlidingDFT_Batch)
    ->Apply(window_power_and_batch_args)
    ->ArgNames({"N", "power_every", "batch"});
BENCHMARK(BM_PocketFFT)
    ->Apply(window_and_power_args)
    ->ArgNames({"N", "power_every"});

int main(int argc, char **argv) {
  if (!check_agreement()) {
    return 1;
  }

  // --csv is ours: strip it before Google Benchmark parses the rest.
  bool csv = false;
  int kept = 1;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--csv") == 0) {
      csv = true;
    } else {
      argv[kept++] = argv[i];
    }
  }
  argc = kept;

  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }
  if (csv) {
    TidyCsvReporter reporter;
    benchmark::RunSpecifiedBenchmarks(&reporter);
  } else {
    benchmark::RunSpecifiedBenchmarks();
  }
  benchmark::Shutdown();
  return 0;
}
