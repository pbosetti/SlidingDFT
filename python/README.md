# Python Wrapper

This folder contains a small `ctypes` wrapper around the native `SlidingDFT`
library. The Python side has no third-party dependencies; it loads a shared
library produced by the project build.

## Build The Native Library

From the repository root:

```sh
cmake -B build
cmake --build build
```

This builds the C bridge library used by Python:

- macOS: `build/libSlidingDFTC.dylib`
- Linux: `build/libSlidingDFTC.so`
- Windows: `build/SlidingDFTC.dll`

The Python wrapper detects automatically whether the native build uses `float`
or `double`.

## Use The Wrapper

From the repository root:

```sh
PYTHONPATH=python python3
```

Then:

```python
from math import pi, sin

from sliding_dft import SlidingDFT

window_size = 256
sdft = SlidingDFT(window_size)

for n in range(5000):
  sample = sin(2.0 * pi * 17.0 * n / window_size)
  sdft.update(sample)

power = sdft.get_power()
print(sdft.bins)
print(power[:8])
```

You can also push a sequence in one call:

```python
sdft.extend([0.0, 1.0, 0.0, -1.0])
```

## Selecting The Shared Library

The wrapper searches these locations by default:

- `build/`
- `build/lib/`
- `python/`
- the exact path from `SLIDINGDFT_LIBRARY`

You can point to a specific build explicitly:

```python
sdft = SlidingDFT(256, library_path="/absolute/path/to/libSlidingDFTC.dylib")
```
