#include <Rcpp.h>

#include <SlidingDFT/sliding_dft.hpp>

#include <span>
#include <type_traits>

namespace {

static_assert(std::is_same_v<data_t, double>,
              "the R package is built with SLIDINGDFT_USE_DOUBLE");

class RSlidingDFT {
public:
  explicit RSlidingDFT(int window_size) : RSlidingDFT(window_size, 4096) {}

  RSlidingDFT(int window_size, int renorm_interval)
      : _sdft(checked_size(window_size, "window_size"),
              checked_size(renorm_interval, "renorm_interval")) {}

  void update(const Rcpp::NumericVector &samples) {
    _sdft.update(std::span<const data_t>(
        samples.begin(), static_cast<std::size_t>(samples.size())));
  }

  Rcpp::NumericVector get_power() const {
    _sdft.get_power(_power);
    return Rcpp::NumericVector(_power.begin(), _power.end());
  }

  int bins() const { return static_cast<int>(_sdft.bins()); }

private:
  static std::size_t checked_size(int value, const char *name) {
    if (value == NA_INTEGER || value <= 0) {
      Rcpp::stop("`%s` must be a positive integer", name);
    }
    return static_cast<std::size_t>(value);
  }

  SlidingDFT _sdft;
  mutable SlidingDFT::vec_t _power;
};

} // namespace

RCPP_MODULE(sliding_dft_module) {
  Rcpp::class_<RSlidingDFT>("SlidingDFT")
      .constructor<int>()
      .constructor<int, int>()
      .method("update", &RSlidingDFT::update,
              "Push a vector of samples into the window")
      .method("get_power", &RSlidingDFT::get_power,
              "Return the current power spectrum")
      .property("bins", &RSlidingDFT::bins, "Number of frequency bins");
}
