# FFT-Based Digital Signal Spectrum Analyzer

> **A from-scratch C++17 implementation of the Cooley-Tukey Radix-2 FFT algorithm for signal frequency analysis**

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://isocpp.org/std/the-standard)
[![CMake](https://img.shields.io/badge/CMake-3.16%2B-green.svg)](https://cmake.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/Tests-18%2F18%20Passed-brightgreen.svg)]()

---

## Project Introduction

This project implements the **Fast Fourier Transform (FFT)** from scratch in C++17, without using any external FFT libraries. It generates a test signal composed of two sinusoids, applies the FFT, and identifies the exact frequencies present in the signal.

The project is designed to be:
- **Interview-ready**: Clean, well-commented code that is easy to explain.
- **Technically correct**: Numerically verified against a direct DFT implementation.
- **Beginner-friendly**: Readable for anyone familiar with basic C++ and data structures.
- **Relevant to DSP/Modem Engineering**: Demonstrates concepts directly used in 4G/5G OFDM modems.

---

## Problem Statement

Given a digital signal sampled at **Fs = 1024 Hz** containing two sinusoids at **50 Hz** and **120 Hz**:

```
x[n] = sin(2π · 50 · n/1024) + 0.5 · sin(2π · 120 · n/1024)
```

**Goal**: Efficiently identify the frequencies present using the FFT algorithm.

**Challenge**: A naive Discrete Fourier Transform (DFT) requires **O(N²)** operations. For N = 1024, that's ~1,000,000 operations. The FFT reduces this to **O(N log N)** ≈ 10,240 operations — roughly **100× faster**.

---

## Technology Stack

| Tool | Purpose |
|------|---------|
| C++17 | Core language |
| STL (`std::complex`, `std::chrono`, `std::vector`) | Data structures, timing |
| CMake 3.16+ | Build system |
| MinGW-w64 GCC | Compiler (Windows) |
| CTest | Automated testing |
| Python + matplotlib | Optional visualization |

---

## How the FFT Algorithm Works

### What is the DFT?

The Discrete Fourier Transform converts a time-domain signal `x[n]` of N samples into N complex frequency-domain values `X[k]`:

```
X[k] = Σ_{n=0}^{N-1}  x[n] · e^(-j·2π·k·n/N)        k = 0, 1, ..., N-1
```

Each `X[k]` tells you how much of frequency `k · Fs / N` Hz is present in the signal. Computing this directly costs **O(N²)**.

### What is the FFT?

The FFT is an algorithm that computes the *exact same result* as the DFT but exploits two mathematical symmetries of the complex exponential to reduce operations to **O(N log₂ N)**:

1. **Periodicity**: `e^(-j·2πk(n+N)/N) = e^(-j·2πkn/N)`
2. **Symmetry**: `e^(-j·2π(k+N/2)n/N) = -e^(-j·2πkn/N)`

This implementation uses the **Cooley-Tukey Radix-2** algorithm, requiring N to be a power of 2.

### Algorithm Steps

```
Step 1: Bit-Reversal Permutation
        Reorder x[n] so index i goes to bit_reverse(i)
        Example (N=8):  [x0,x1,x2,x3,x4,x5,x6,x7]
                     → [x0,x4,x2,x6,x1,x5,x3,x7]

Step 2: log₂(N) Butterfly Passes
        For each stage s = 1, 2, ..., log₂N:
          For each group of 2^s elements:
            For each pair (top, bottom):
              t = W · bottom
              top    = top + t    ← butterfly upper output
              bottom = top - t    ← butterfly lower output
              W = W · W_m         ← advance twiddle factor

Result: x[k] now contains X[k] (FFT output)
```

### Why Bit-Reversal?

The recursive Cooley-Tukey algorithm repeatedly splits the input into **even** and **odd** indexed elements. When you unroll all the recursive splitting, the data lands in bit-reversed order. By pre-permuting the input before the butterfly passes, we can run all butterflies in-place in a simple loop.

### Why O(N log N)?

- There are **log₂(N)** butterfly passes
- Each pass processes **N/2** butterfly pairs
- Each butterfly = 1 complex multiply + 2 complex adds = **O(1)**
- Total: **(N/2) × log₂(N)** operations = **O(N log N)**

### DFT vs FFT Complexity Table

| N | DFT (N²) | FFT (N log₂N) | Speedup |
|---|----------|----------------|---------|
| 64 | 4,096 | 384 | ~10.7× |
| 256 | 65,536 | 2,048 | ~32× |
| 1,024 | 1,048,576 | 10,240 | ~102× |
| 4,096 | 16,777,216 | 49,152 | ~341× |
| 65,536 | 4,294,967,296 | 1,048,576 | ~4,096× |

---

## Architecture

```
FFT-Spectrum-Analyzer/
│
├── include/
│   ├── fft.h               ← FFT/IFFT/DFT declarations + full DSP documentation
│   └── signal_analyzer.h   ← Signal generation and analysis declarations
│
├── src/
│   ├── fft.cpp             ← Iterative Cooley-Tukey FFT (core algorithm)
│   ├── signal_analyzer.cpp ← Signal generator, spectrum analysis, CSV export
│   └── main.cpp            ← Entry point, CLI, timing, display
│
├── tests/
│   └── test_fft.cpp        ← 18-test suite (impulse, DC, sinusoid, IFFT, invalid)
│
├── scripts/
│   └── plot.py             ← matplotlib visualization (optional)
│
├── results/                ← Generated CSV files land here
│   ├── signal.csv
│   ├── spectrum.csv
│   └── spectrum_plot.png
│
└── CMakeLists.txt          ← Build system (two targets + CTest)
```

### Data Flow

```
[main.cpp]
    │
    ├─ parseArgs()           → SignalParams {Fs, N, f1, f2, amp1, amp2}
    │
    ├─ validateParams()      → Checks Nyquist, power-of-2, positivity
    │
    ├─ generateSignal()      → ComplexVec [N samples]
    │      x[n] = A1·sin(2πf1·n/Fs) + A2·sin(2πf2·n/Fs)
    │
    ├─ fft::fft()            → In-place FFT → ComplexVec X[k]
    │      ├─ bitReversalPermutation()
    │      └─ log₂N butterfly passes
    │
    ├─ computeMagnitudeSpectrum()  → vector<double> magnitudes[N/2+1]
    │
    ├─ detectPeaks()         → vector<FrequencyPeak> {freq, amp, bin}
    │
    └─ exportSignalToCSV() + exportSpectrumToCSV()  → results/*.csv
```

---

## Build and Execution (Windows PowerShell)

### Prerequisites

1. **CMake** ≥ 3.16: [cmake.org/download](https://cmake.org/download/)
2. **MinGW-w64 GCC** (via MSYS2): [msys2.org](https://www.msys2.org/)
   ```powershell
   # Install MSYS2, then in MSYS2 terminal:
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake
   ```
3. Add `C:\msys64\mingw64\bin` to your Windows PATH.

### Build

```powershell
# Navigate to the project directory
cd FFT-Spectrum-Analyzer

# Configure (creates build/ directory with MinGW Makefiles)
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release

# Build all targets (spectrum_analyzer.exe + test_fft.exe)
cmake --build build --config Release
```

### Run the Analyzer

```powershell
# Default: Fs=1024, N=1024, f1=50Hz, f2=120Hz
.\build\spectrum_analyzer.exe

# Custom: Fs=2048, N=2048, f1=100Hz, f2=300Hz, amp1=1.0, amp2=0.8
.\build\spectrum_analyzer.exe 2048 2048 100 300 1.0 0.8
```

### Run Tests

```powershell
ctest --test-dir build --verbose
```

### Visualize Results (Optional)

```powershell
pip install matplotlib numpy
python scripts\plot.py
```

---

## Example Terminal Output

```
╔══════════════════════════════════════════════╗
║   FFT-Based Digital Signal Spectrum Analyzer  ║
║   Built with C++17 | Cooley-Tukey Radix-2    ║
╚══════════════════════════════════════════════╝

Configuration
─────────────────────────────────────
  Sampling Frequency (Fs) : 1024.0 Hz
  Number of Samples  (N)  : 1024
  Signal Frequency 1 (f1) : 50.0 Hz  (amplitude 1.0)
  Signal Frequency 2 (f2) : 120.0 Hz  (amplitude 0.5)
  Frequency Resolution    : 1.0 Hz/bin
  Nyquist Frequency       : 512.0 Hz
─────────────────────────────────────

Generating signal...
  Generated 1024 samples.

Performing FFT...
  FFT execution time: 158 µs
  Theoretical ops: FFT O(N log N) ≈ 10240 | DFT O(N²) ≈ 1048576

Detected Frequencies
─────────────────────────────────────
  Peak 1: 50.0 Hz  (normalized amplitude = 1.0000,  bin index = 50)
  Peak 2: 120.0 Hz  (normalized amplitude = 0.5000,  bin index = 120)
─────────────────────────────────────

[DFT vs FFT Comparison on first 64 samples]
  DFT time  : 558 µs  (O(N²), N=64)
  FFT time  : 7 µs  (O(N log N), N=64)
  Max |DFT[k] - FFT[k]| = 4.09e-13  (near machine epsilon ~1e-13)
  ✓ FFT output matches DFT output

Exporting results...
  Signal   → results/signal.csv
  Spectrum → results/spectrum.csv

╔══════════════════════════════════════════════╗
║               FFT COMPLETE ✓                 ║
╚══════════════════════════════════════════════╝
```

---

## Generated Graphs (plot.py)

The Python script produces two subplots on a dark-themed canvas:

**Plot 1 — Time-Domain Signal**
Shows the first 100 ms of the composite waveform `x(t) = sin(2π·50·t) + 0.5·sin(2π·120·t)`. You can see the characteristic beat pattern where the two sinusoids reinforce and partially cancel each other.

**Plot 2 — One-Sided FFT Magnitude Spectrum**
Shows normalized amplitude vs. frequency (0 to 512 Hz). Two sharp vertical spikes appear at exactly **50 Hz** and **120 Hz**, annotated with arrows. The rest of the spectrum is near zero (no noise in this simulated signal).

---

## Testing

```powershell
ctest --test-dir build --verbose
```

| # | Test | What it Verifies |
|---|------|-----------------|
| 1 | isPowerOfTwo utility | Correct power-of-2 detection |
| 2 | Impulse FFT | X[k] = 1 for all k (N=8) |
| 3 | DC signal FFT | X[0] = N·C, X[k≠0] = 0 |
| 4 | 50 Hz sinusoid | Peak at bin 50, magnitude ≈ N/2 |
| 5 | IFFT reconstruction | max|IFFT(FFT(x)) - x| < 1e-9 |
| 6 | Invalid input | Exception for N=100 (non-power of 2) |
| 7 | Full pipeline (50+120 Hz) | Both peaks detected correctly |

**Result: 18/18 tests passed.**

---

## Project Limitations

1. **Bin-aligned frequencies only**: Peak detection works perfectly when `f × N / Fs` is an integer (a "bin-aligned" frequency). For non-bin-aligned frequencies, **spectral leakage** spreads energy across adjacent bins, causing amplitude estimation errors. Real-world systems use **windowing functions** (Hanning, Blackman) to reduce leakage.

2. **No noise**: This project simulates a noiseless signal. With additive noise, naive peak detection can fail. Production systems use **CFAR (Constant False Alarm Rate)** or more robust algorithms.

3. **No real-time input**: This is an offline batch analyzer. Embedded/real-time systems use DMA-based ADC pipelines and process signal blocks in a streaming fashion.

4. **Fixed radix-2**: Only supports N = power of 2. Mixed-radix FFT (e.g., FFTW) handles arbitrary N.

5. **Not embedded firmware**: This runs on a PC with an OS. Embedded FFT for microcontrollers requires fixed-point arithmetic, SIMD intrinsics, and tight memory budgets.

---

## Relevance to Cellular Modem PHY Firmware

This project directly demonstrates skills relevant to the Apple Cellular Modem PHY team:

| This Project | Modem PHY Application |
|---|---|
| FFT frequency detection | OFDM symbol demodulation (4G LTE, 5G NR) |
| Iterative Cooley-Tukey FFT | Basis of IFFT at OFDM transmitter |
| Bin-to-Hz conversion | Subcarrier frequency mapping |
| `std::complex<double>` | Complex baseband signal processing |
| Amplitude normalization | Channel equalization after FFT |
| IFFT reconstruction test | Transmitter-receiver loopback validation |
| Numerical correctness with tolerance | Floating-point guard bands in hardware/software co-design |

In LTE, the eNodeB transmitter uses a **2048-point IFFT** (at 20 MHz bandwidth) to convert frequency-domain OFDM symbols to a time-domain waveform. The UE receiver uses a **2048-point FFT** to demodulate. This runs thousands of times per second and must be implemented in highly optimized C/C++ or RTL on the modem DSP core.

---

## License

MIT License — free to use, modify, and distribute.
