#include <SlidingDFT/c_api.h>
#include <SlidingDFT/sliding_dft.hpp>

#include <algorithm>
#include <exception>
#include <memory>
#include <span>
#include <string>

struct SlidingDFTHandle {
  explicit SlidingDFTHandle(std::size_t window_size, std::size_t renorm_interval)
      : instance(window_size, renorm_interval) {}

  SlidingDFT instance;
  // Reused by sliding_dft_get_power() to avoid an allocation per call.
  mutable SlidingDFT::vec_t power;
};

namespace {

thread_local std::string LastError;

void set_error(const char *message) { LastError = message; }

void set_error(const std::exception &error) { LastError = error.what(); }

} // namespace

extern "C" {

SlidingDFTHandle *sliding_dft_create(std::size_t window_size,
                                     std::size_t renorm_interval) {
  try {
    LastError.clear();
    return new SlidingDFTHandle(window_size, renorm_interval);
  } catch (const std::exception &error) {
    set_error(error);
    return nullptr;
  } catch (...) {
    set_error("unknown error");
    return nullptr;
  }
}

void sliding_dft_destroy(SlidingDFTHandle *handle) { delete handle; }

int sliding_dft_update(SlidingDFTHandle *handle, data_t sample) {
  if (handle == nullptr) {
    set_error("null SlidingDFT handle");
    return -1;
  }

  try {
    LastError.clear();
    handle->instance.update(sample);
    return 0;
  } catch (const std::exception &error) {
    set_error(error);
    return -1;
  } catch (...) {
    set_error("unknown error");
    return -1;
  }
}

int sliding_dft_update_many(SlidingDFTHandle *handle, const data_t *samples,
                            std::size_t count) {
  if (handle == nullptr) {
    set_error("null SlidingDFT handle");
    return -1;
  }
  if (samples == nullptr && count > 0) {
    set_error("null sample buffer");
    return -1;
  }

  try {
    LastError.clear();
    handle->instance.update(std::span(samples, count));
    return 0;
  } catch (const std::exception &error) {
    set_error(error);
    return -1;
  } catch (...) {
    set_error("unknown error");
    return -1;
  }
}

std::size_t sliding_dft_bins(const SlidingDFTHandle *handle) {
  if (handle == nullptr) {
    set_error("null SlidingDFT handle");
    return 0;
  }

  LastError.clear();
  return handle->instance.bins();
}

std::size_t sliding_dft_sample_size() { return sizeof(data_t); }

int sliding_dft_get_power(const SlidingDFTHandle *handle, data_t *out,
                          std::size_t out_size) {
  if (handle == nullptr) {
    set_error("null SlidingDFT handle");
    return -1;
  }
  if (out == nullptr) {
    set_error("null output buffer");
    return -1;
  }

  try {
    LastError.clear();

    SlidingDFT::vec_t &power = handle->power;
    handle->instance.get_power(power);
    if (out_size < power.size()) {
      set_error("output buffer too small");
      return -1;
    }

    std::copy(power.begin(), power.end(), out);
    return 0;
  } catch (const std::exception &error) {
    set_error(error);
    return -1;
  } catch (...) {
    set_error("unknown error");
    return -1;
  }
}

const char *sliding_dft_last_error() { return LastError.c_str(); }

}
