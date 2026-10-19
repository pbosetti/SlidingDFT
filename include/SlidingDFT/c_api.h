#pragma once

#include <SlidingDFT/types.hpp>

#include <cstddef>

#if defined(_WIN32)
#if defined(SLIDINGDFT_C_EXPORTS)
#define SLIDINGDFT_C_API __declspec(dllexport)
#else
#define SLIDINGDFT_C_API __declspec(dllimport)
#endif
#else
#define SLIDINGDFT_C_API
#endif

extern "C" {

/**
 * @struct SlidingDFTHandle
 * @brief Opaque handle for the C interface to `SlidingDFT`.
 *
 * Instances are created with sliding_dft_create() and released with
 * sliding_dft_destroy().
 */
struct SlidingDFTHandle;

/**
 * @brief Create a new sliding DFT instance.
 * @param window_size Number of samples in the analysis window.
 * @param renorm_interval Number of updates between exact recomputations of
 *        the bins.
 * @return A valid handle on success, or `nullptr` on error.
 *
 * When creation fails, call sliding_dft_last_error() to retrieve the reason.
 */
SLIDINGDFT_C_API SlidingDFTHandle *sliding_dft_create(std::size_t window_size,
                                                      std::size_t renorm_interval);

/**
 * @brief Destroy a sliding DFT instance created by sliding_dft_create().
 * @param handle Handle to destroy. `nullptr` is allowed.
 */
SLIDINGDFT_C_API void sliding_dft_destroy(SlidingDFTHandle *handle);

/**
 * @brief Push one new sample into the analysis window.
 * @param handle Sliding DFT handle.
 * @param sample New input sample.
 * @return `0` on success, `-1` on error.
 *
 * On failure, call sliding_dft_last_error() for details.
 */
SLIDINGDFT_C_API int sliding_dft_update(SlidingDFTHandle *handle, data_t sample);

/**
 * @brief Push a block of samples into the analysis window, oldest first.
 * @param handle Sliding DFT handle.
 * @param samples Input samples, `count` elements.
 * @param count Number of samples; `0` is allowed.
 * @return `0` on success, `-1` on error.
 *
 * Gives the same result as calling sliding_dft_update() on each sample, but is
 * much faster. On failure, call sliding_dft_last_error() for details.
 */
SLIDINGDFT_C_API int sliding_dft_update_many(SlidingDFTHandle *handle,
                                             const data_t *samples,
                                             std::size_t count);

/**
 * @brief Return the number of frequency bins for an instance.
 * @param handle Sliding DFT handle.
 * @return The number of bins, or `0` on error.
 *
 * The bin count is `window_size / 2 + 1`.
 */
SLIDINGDFT_C_API std::size_t sliding_dft_bins(const SlidingDFTHandle *handle);

/**
 * @brief Return the size in bytes of the sample type used by the library.
 * @return `sizeof(data_t)` for the current native build.
 *
 * This allows foreign-language bindings to determine whether the library was
 * built with `float` or `double` samples.
 */
SLIDINGDFT_C_API std::size_t sliding_dft_sample_size();

/**
 * @brief Compute the current power spectrum.
 * @param handle Sliding DFT handle.
 * @param out Output buffer receiving `sliding_dft_bins(handle)` elements.
 * @param out_size Number of elements available in @p out.
 * @return `0` on success, `-1` on error.
 *
 * The caller owns the output buffer. On failure, call
 * sliding_dft_last_error() for details.
 */
SLIDINGDFT_C_API int sliding_dft_get_power(const SlidingDFTHandle *handle,
                                           data_t *out,
                                           std::size_t out_size);

/**
 * @brief Return the last error message for the current thread.
 * @return A pointer to a null-terminated error string.
 *
 * The returned pointer remains owned by the library and must not be freed by
 * the caller. Its contents may change after the next API call on the same
 * thread.
 */
SLIDINGDFT_C_API const char *sliding_dft_last_error();

}
