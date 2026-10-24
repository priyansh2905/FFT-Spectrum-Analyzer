# INTERVIEW_PREP.md — FFT-Based Digital Signal Spectrum Analyzer

> **Apple Firmware Development Engineering Intern — Cellular Modem PHY Team**
> Interview preparation guide in plain English.

---

## PART A: Project Explanation

---

### 30-Second Elevator Pitch

> "I built an FFT-Based Digital Signal Spectrum Analyzer in C++17 from scratch. The program generates a composite signal containing two frequencies, then uses the Cooley-Tukey Fast Fourier Transform — which I implemented myself without any library — to analyze the signal and identify those frequencies. The FFT runs in O(N log N) time compared to O(N²) for a naive approach. I also wrote a test suite with 18 tests and a Python script to visualize the results. The core algorithm is directly used in OFDM modems like 4G LTE."

---

### 2-Minute Interview Introduction

> "For my interview project, I implemented an FFT-based spectrum analyzer in C++17.
>
> The project has three main parts: signal generation, FFT computation, and frequency analysis.
>
> First, I generate a time-domain signal containing two sinusoids — 50 Hz and 120 Hz — sampled at 1024 Hz. This models a simplified radio signal.
>
> Second, I apply the Fast Fourier Transform. I implemented the iterative Cooley-Tukey Radix-2 FFT from scratch using bit-reversal permutation and butterfly operations. This transforms the signal from the time domain to the frequency domain in O(N log N) time.
>
> Third, I analyze the FFT output, compute the magnitude spectrum, normalize it correctly for a one-sided display, and identify the dominant frequency peaks.
>
> The program correctly detects 50 Hz and 120 Hz. I verified correctness by comparing the FFT output against a direct DFT computation — they match to within 4 × 10⁻¹³, which is essentially machine precision.
>
> I also integrated CTest for automated testing with 18 tests covering impulse response, DC signal, sinusoid detection, inverse FFT reconstruction, and error handling.
>
> The reason I built this is that FFT is the cornerstone of OFDM — the modulation scheme used in 4G LTE and 5G NR. Understanding FFT deeply is essential for PHY-layer firmware work."

---

### Step-by-Step Processing Pipeline

```
Step 1: Parse CLI Arguments
        User provides: Fs, N, f1, f2, amp1, amp2
        Default: Fs=1024, N=1024, f1=50, f2=120, amp1=1.0, amp2=0.5

Step 2: Validate Parameters
        • N must be a power of 2
        • Fs > 2 × max(f1, f2)   [Nyquist criterion]
        • f1, f2, amp1, amp2 > 0

Step 3: Generate Signal
        For each sample n = 0 to N-1:
          x[n] = amp1 × sin(2π × f1 × n/Fs) + amp2 × sin(2π × f2 × n/Fs)
        Store as complex<double> (imaginary part = 0)

Step 4: FFT (In-Place Cooley-Tukey)
        4a. Bit-reversal permutation: reorder x[n] by bit-reversed indices
        4b. log₂(N) butterfly passes:
            Each pass doubles the sub-DFT size (1→2→4→...→N)
            Each butterfly: (top, bottom) → (top + W×bottom, top - W×bottom)
        After all passes: x[k] = X[k] (FFT output)

Step 5: Magnitude Spectrum
        For each bin k = 0 to N/2:
          magnitude[k] = |X[k]| (normalized for one-sided spectrum)

Step 6: Peak Detection
        Find all local maxima above threshold
        Sort by amplitude, return top 2
        Convert bin index → Hz: freq = k × Fs / N

Step 7: Display + Export
        Print detected frequencies to console
        Write signal.csv and spectrum.csv to results/
```

---

### Source File Guide

#### `include/fft.h`
Declares the FFT public API: `fft()`, `ifft()`, `dft()`, `isPowerOfTwo()`, `bitReversalPermutation()`. Contains extensive documentation explaining DFT, FFT, butterfly operations, bit-reversal, complexity, and OFDM relevance.

#### `include/signal_analyzer.h`
Declares the signal processing pipeline: `SignalParams` struct, `FrequencyPeak` struct, `generateSignal()`, `computeMagnitudeSpectrum()`, `detectPeaks()`, `exportSignalToCSV()`, `exportSpectrumToCSV()`, `validateParams()`.

#### `src/fft.cpp`
The heart of the project. Implements:
- `isPowerOfTwo(n)`: Bitwise check `n & (n-1) == 0`
- `bitReversalPermutation(x)`: Iterative bit reversal by bit-by-bit construction
- `fft(x)`: Three nested loops — groups, butterfly pairs, twiddle advance
- `ifft(x)`: Conjugation trick — conj(FFT(conj(X))) / N
- `dft(x)`: Naive O(N²) for comparison

#### `src/signal_analyzer.cpp`
Implements `generateSignal()` using the sinusoid formula, `computeMagnitudeSpectrum()` with proper one-sided normalization, `detectPeaks()` using local-maxima search, and two CSV export functions.

#### `src/main.cpp`
Orchestrates the entire pipeline: parses arguments, validates, generates, runs FFT, measures time with `std::chrono`, displays results, runs DFT comparison, exports CSVs.

#### `tests/test_fft.cpp`
Self-contained test runner (no external framework). Defines a `check(condition, name)` helper and seven test functions. Returns exit code 0 if all pass (CTest requirement).

---

## PART B: Core Technical Concepts

---

### 1. What is DSP?

**DSP (Digital Signal Processing)** means processing signals (audio, radio, video) using numbers on a computer, instead of analog circuits.

**Example**: Your phone's microphone captures sound waves. These waves are converted to numbers by an ADC (Analog-to-Digital Converter). DSP algorithms then process these numbers — to cancel noise, compress audio, or decode a phone call.

**Key idea**: Work with signals as arrays of numbers.

---

### 2. What is DFT?

**DFT (Discrete Fourier Transform)** converts a time-domain signal into a frequency-domain representation.

Think of it like this: you hear a chord played on a piano. You can't easily tell which individual notes are being played just by looking at the sound wave. The DFT "decomposes" the chord into individual notes (frequencies).

**Formula**:
```
X[k] = Σ_{n=0}^{N-1}  x[n] · e^(-j·2π·k·n/N)
```

- `x[n]` = time-domain sample at index n
- `X[k]` = frequency-domain value at bin k
- Output `X[k]` is complex: contains both amplitude and phase at frequency `k · Fs/N` Hz

**Complexity**: O(N²) — slow for large N.

---

### 3. What is FFT?

**FFT (Fast Fourier Transform)** is an *algorithm* that computes the *exact same result* as DFT but much faster by exploiting mathematical symmetry.

**Analogy**: Instead of computing all N² combinations, the FFT "reuses" previously computed results. It's like solving a complex puzzle by solving smaller identical sub-puzzles and combining them.

**Complexity**: O(N log₂N) — much faster than O(N²).

---

### 4. Why is FFT Faster?

The DFT of N points can be split into:
- DFT of the **even-indexed** samples (N/2 points)
- DFT of the **odd-indexed** samples (N/2 points)
- Then combine them in N/2 "butterfly" operations

This is the divide-and-conquer idea. Each half can be split again, until you reach 1-point DFTs (which are trivial). After log₂N levels of splitting, you only need (N/2) × log₂N operations total.

**Numbers for N = 1024**:
- DFT: 1,048,576 operations
- FFT: 10,240 operations
- **~100× faster**

---

### 5. What is a Butterfly Operation?

The butterfly is the fundamental building block of FFT.

It takes two complex numbers `a` and `b`, plus a "twiddle factor" `W`:

```
a' = a + W × b
b' = a - W × b
```

**Why "butterfly"?** Draw the signal flow diagram — lines cross from `a` and `b` to `a'` and `b'`, forming a shape that looks like a butterfly.

**Small example (N=2, W=1)**:
```
a = 3,  b = 1,  W = 1
a' = 3 + 1×1 = 4
b' = 3 - 1×1 = 2
This is the 2-point DFT!
```

---

### 6. What is Bit Reversal?

When we split the DFT recursively (even/odd, even/odd...), the final ordering of inputs is the bit-reversed order of their original indices.

**Example for N=8 (3 bits)**:

| Original Index | Binary | Bit-Reversed Binary | Bit-Reversed Index |
|---|---|---|---|
| 0 | 000 | 000 | 0 |
| 1 | 001 | 100 | 4 |
| 2 | 010 | 010 | 2 |
| 3 | 011 | 110 | 6 |
| 4 | 100 | 001 | 1 |
| 5 | 101 | 101 | 5 |
| 6 | 110 | 011 | 3 |
| 7 | 111 | 111 | 7 |

**Why?** Without bit-reversal, the butterfly outputs go to the wrong positions. Bit-reversal permutation "puts everything in the right place" before the butterfly passes begin, enabling an in-place, iterative computation.

---

### 7. What is Sampling Frequency?

**Sampling frequency (Fs)** is how many times per second an analog signal is measured (sampled) to convert it to digital.

**Example**: Fs = 1024 Hz means we take 1024 measurements per second. Between each measurement, `T = 1/Fs = 0.976 ms` passes.

**Higher Fs** → Can represent higher frequencies → Better signal quality → More data.

---

### 8. What is the Nyquist Theorem?

**Nyquist Theorem**: To correctly represent a signal containing frequency `f`, you must sample at least at `Fs = 2f`.

More precisely, to avoid aliasing: **Fs > 2 × f_max**.

**In this project**: f_max = 120 Hz, Fs = 1024 Hz. Since 1024 > 2 × 120 = 240, we satisfy Nyquist. The Nyquist frequency (highest representable) = Fs/2 = 512 Hz.

---

### 9. What is Aliasing?

**Aliasing** happens when you sample too slowly (Fs ≤ 2f). High-frequency components "fold back" and appear as lower frequencies — a false (aliased) frequency appears in your spectrum.

**Example**: If Fs = 100 Hz and your signal has a 90 Hz component, it appears as 100 - 90 = 10 Hz in the FFT output. This is wrong — it's an alias.

**Prevention**: Use an anti-aliasing lowpass filter (LPF) before the ADC to remove all frequencies above Fs/2.

---

### 10. What is a Frequency Bin?

The FFT outputs N complex values `X[0], X[1], ..., X[N-1]`. Each output `X[k]` corresponds to a specific frequency:

```
frequency of bin k = k × Fs / N   (in Hz)
```

**Example** (N=1024, Fs=1024):
- Bin 0 → 0 Hz (DC)
- Bin 1 → 1 Hz
- Bin 50 → 50 Hz ← our first sinusoid
- Bin 120 → 120 Hz ← our second sinusoid
- Bin 512 → 512 Hz (Nyquist, last useful bin)

**Frequency resolution** = Fs/N = 1024/1024 = 1 Hz/bin. You can only detect frequencies at multiples of 1 Hz with this setup.

---

### 11. Time Domain vs Frequency Domain

| Aspect | Time Domain | Frequency Domain |
|---|---|---|
| Question answered | What is the signal at each moment? | What frequencies are in the signal? |
| x-axis | Time (seconds) | Frequency (Hz) |
| y-axis | Amplitude | Magnitude (strength of each frequency) |
| Representation | Waveform | Spectrum |
| Useful for | Understanding signal shape | Understanding signal composition |

**Analogy**: Time domain = the sound of a chord. Frequency domain = the piano sheet music (which individual notes).

---

### 12. Why Complex Numbers in FFT?

The DFT formula uses `e^(-j·2π·k·n/N)` — a complex exponential. By Euler's formula:

```
e^(-jθ) = cos(θ) - j·sin(θ)
```

This simultaneously tracks both **cosine** and **sine** components of each frequency, capturing **both amplitude and phase** in a single complex number.

Without complex numbers, you'd need two separate calculations for sine and cosine components. Using `std::complex<double>` makes the butterfly formula elegant and matches the mathematical definition exactly.

For a real-valued input signal, `X[k]` and `X[N-k]` are complex conjugates — that's why we only need to look at the first N/2 bins (the "one-sided spectrum").

---

### 13. Why is FFT Useful in OFDM?

**OFDM (Orthogonal Frequency Division Multiplexing)** is used in 4G LTE, 5G NR, WiFi, and DVB.

**Concept**: Instead of sending one wide-band signal, split data into many narrow sub-channels (subcarriers) at different frequencies. Each subcarrier is orthogonal to all others (no interference).

**How FFT fits in**:
```
Transmitter:
  Data bits → assign to subcarriers (frequency domain) → IFFT → time-domain waveform → RF

Receiver:
  RF → ADC → FFT → frequency domain → read data from each subcarrier → decode bits
```

The IFFT/FFT pair is what makes OFDM practical. Without FFT's O(N log N) efficiency, running a 2048-point FFT thousands of times per second in a modem chip would be impossible.

**In 4G LTE**: 2048-point FFT at 30.72 MHz sampling rate, ~14,000 times per second.

---

### 14. Time and Space Complexity

| Algorithm | Time Complexity | Space Complexity |
|---|---|---|
| Direct DFT | O(N²) | O(N) |
| Cooley-Tukey FFT | O(N log₂ N) | O(N) — in-place |
| IFFT (conjugation method) | O(N log₂ N) | O(N) — in-place |

**Space**: The in-place FFT modifies the input array directly, requiring no extra memory beyond the input buffer. This is critical for embedded systems with limited RAM.

---

### 15. Adapting for an Embedded Microcontroller

1. **Fixed-point arithmetic**: Replace `double` with `int16_t` or `int32_t`. Use Q-format (e.g., Q15) arithmetic. Reduces memory and eliminates FPU dependency.

2. **Pre-computed twiddle table**: Compute twiddle factors at compile time (stored in flash). Avoids `sin()` and `cos()` at runtime.

3. **DMA-driven ADC**: Use DMA to fill the sample buffer while the CPU processes the previous block — pipelining.

4. **SIMD/DSP intrinsics**: Use ARM CMSIS-DSP library (which provides `arm_cfft_f32()`) for hardware-accelerated butterfly on Cortex-M4/M7 with FPU.

5. **Cache alignment**: Align the sample buffer to cache line boundaries to maximize memory access efficiency.

6. **Static allocation**: Replace `std::vector` with fixed-size arrays to avoid heap allocation in ISRs (interrupt service routines).

---

## PART C: Coding Walkthrough

---

### 4-Point FFT Dry Run

**Input**: x = [1, 2, 3, 4]  (N=4, log₂4 = 2 stages)

**Step 1: Bit-Reversal (2 bits)**

| Index | Binary | Reversed | Reversed Index |
|---|---|---|---|
| 0 | 00 | 00 | 0 |
| 1 | 01 | 10 | 2 |
| 2 | 10 | 01 | 1 |
| 3 | 11 | 11 | 3 |

Swap x[1] ↔ x[2]:  `x = [1, 3, 2, 4]`

---

**Step 2: Stage 1 (m=2, half_m=1)**

Twiddle base: `W₂ = e^(-j·2π/2) = e^(-jπ) = -1`

**Group 0 (indices 0,1)**: j=0, W = 1
- t = 1 × x[1] = 1 × 3 = 3
- x[1] = x[0] - t = 1 - 3 = **-2**
- x[0] = x[0] + t = 1 + 3 = **4**

**Group 1 (indices 2,3)**: j=0, W = 1
- t = 1 × x[3] = 1 × 4 = 4
- x[3] = x[2] - t = 2 - 4 = **-2**
- x[2] = x[2] + t = 2 + 4 = **6**

After Stage 1: `x = [4, -2, 6, -2]`

---

**Step 3: Stage 2 (m=4, half_m=2)**

Twiddle base: `W₄ = e^(-j·2π/4) = e^(-jπ/2) = cos(-π/2) + j·sin(-π/2) = 0 - j`

**Group 0 (indices 0,1,2,3):**

j=0: W = 1
- t = 1 × x[2] = 1 × 6 = 6
- x[2] = x[0] - t = 4 - 6 = **-2**
- x[0] = x[0] + t = 4 + 6 = **10**

j=1: W = W × W₄ = 1 × (0-j) = **-j**
- t = (-j) × x[3] = (-j) × (-2) = 2j
- x[3] = x[1] - t = -2 - 2j = **-2 - 2j**
- x[1] = x[1] + t = -2 + 2j = **-2 + 2j**

---

**Final FFT Result**: `X = [10, -2+2j, -2, -2-2j]`

**Verification** (manual DFT for k=0):
- X[0] = x[0]+x[1]+x[2]+x[3] = 1+2+3+4 = **10** ✓

**Verification** (k=1):
- X[1] = 1 + 2·e^(-jπ/2) + 3·e^(-jπ) + 4·e^(-j3π/2)
       = 1 + 2·(-j) + 3·(-1) + 4·(j)
       = 1 - 2j - 3 + 4j = **-2 + 2j** ✓

---

### C++ Code Walkthrough — `fft()` Function

```cpp
void fft(ComplexVec& x) {
```
> Takes the input array by reference — will modify it in-place.

```cpp
    const std::size_t N = x.size();
    if (!isPowerOfTwo(N)) {
        throw std::invalid_argument("...");
    }
```
> Get size, validate. Throws immediately if N isn't a power of 2.

```cpp
    bitReversalPermutation(x);
```
> Reorder array by bit-reversed indices. After this, the data is in the order that butterflies expect.

```cpp
    for (std::size_t s = 1; (1u << s) <= N; ++s) {
```
> Loop over stages: s = 1, 2, ..., log₂N.
> `1u << s` = 2^s. The loop runs while 2^s ≤ N, i.e., log₂N times.

```cpp
        const std::size_t m = 1u << s;        // m = 2^s (sub-DFT size)
        const std::size_t half_m = m >> 1;    // m/2 = 2^(s-1)
```
> `m` = size of DFT being combined this stage.
> `half_m` = number of butterfly pairs per group.

```cpp
        const Complex W_m = Complex(
            std::cos(-2.0 * PI / m),
            std::sin(-2.0 * PI / m)
        );
```
> `W_m = e^(-j·2π/m)` — the principal m-th root of unity. This is the twiddle factor "step" for this stage.

```cpp
        for (std::size_t k = 0; k < N; k += m) {
```
> Outer loop: process each group of m elements.
> Groups start at k=0, m, 2m, ... until all N elements are covered.

```cpp
            Complex W = Complex(1.0, 0.0);
```
> Start with W = W_m^0 = 1.

```cpp
            for (std::size_t j = 0; j < half_m; ++j) {
```
> Inner loop: j-th butterfly pair within this group.

```cpp
                const Complex t = W * x[idx_bot];
                x[idx_bot] = x[idx_top] - t;
                x[idx_top] = x[idx_top] + t;
```
> The butterfly! Computes the two outputs using the twiddle-weighted "bottom" element.

```cpp
                W *= W_m;
```
> Advance twiddle: W_m^j × W_m = W_m^(j+1). Avoids calling cos/sin again.

---

### `bitReversalPermutation()` Walkthrough

```cpp
std::size_t bits = 0;
std::size_t temp = N;
while (temp > 1) { temp >>= 1; ++bits; }
```
> Count how many bits are needed: for N=8, bits=3.

```cpp
for (std::size_t i = 0; i < N; ++i) {
    std::size_t j = 0;
    std::size_t val = i;
    for (std::size_t b = 0; b < bits; ++b) {
        j = (j << 1) | (val & 1);   // Append LSB of val to j
        val >>= 1;                    // Shift val right
    }
    if (j > i) std::swap(x[i], x[j]);
}
```
> For each index i, build its bit-reversed index j one bit at a time.
> Only swap when j > i to avoid swapping the same pair twice.

---

## PART D: 20 Interview Questions & Answers

---

### Q1: What is the time complexity of your FFT and why?

**A**: O(N log₂ N).

The algorithm has log₂N stages. In each stage, we perform N/2 butterfly operations. Each butterfly is O(1) — one complex multiply and two complex adds. Total: (N/2) × log₂N = O(N log N).

---

### Q2: Why does your FFT require N to be a power of 2?

**A**: The Cooley-Tukey Radix-2 algorithm splits N-point DFT into two N/2-point DFTs. For this to work all the way down to 1-point DFTs, N must be divisible by 2 at every level, which requires N = 2^k.

Mixed-radix FFTs (like those in FFTW) handle arbitrary N by splitting into different-sized prime factors, but they are more complex to implement.

---

### Q3: What is the purpose of bit-reversal permutation?

**A**: The recursive Cooley-Tukey algorithm splits the input into even and odd indices repeatedly. When you fully unroll all recursion levels, the natural ordering of the data becomes bit-reversed. Bit-reversal permutation pre-arranges the data into this order before the butterfly passes, enabling an efficient in-place iterative computation without recursion overhead.

---

### Q4: What is a twiddle factor and where does it come from?

**A**: A twiddle factor is a complex exponential `W_N^k = e^(-j·2πk/N)`. It comes from the DFT formula `e^(-j·2πkn/N)`. When we combine two N/2-point DFTs using the Cooley-Tukey butterfly, the factor `e^(-j·2πk/N)` that scales the odd DFT is the twiddle factor. It "rotates" the odd-DFT output by the correct phase angle before combining.

---

### Q5: How does your IFFT work?

**A**: I use the conjugation method:

```
IFFT(X) = conj(FFT(conj(X))) / N
```

Taking the conjugate of the input flips the sign of the exponent in the DFT formula from `-j` to `+j`, which is the IFFT formula. After applying forward FFT, conjugate again and divide by N for normalization. This reuses the FFT code without a separate IFFT implementation.

---

### Q6: How do you convert an FFT bin index to a frequency in Hz?

**A**: `frequency_Hz = k × Fs / N`

Where k = bin index, Fs = sampling frequency, N = number of samples.

**Example**: k=50, Fs=1024, N=1024 → frequency = 50 × 1024/1024 = **50 Hz**.

**Intuition**: The FFT analyzes one complete block of duration `T = N/Fs` seconds. The minimum frequency that completes exactly one cycle in time T is `1/T = Fs/N` Hz. Bin k represents k cycles in time T.

---

### Q7: Why do you only look at N/2 bins for a real-valued signal?

**A**: For a real-valued input signal, the FFT output is **conjugate symmetric**: `X[k] = conj(X[N-k])`. This means bins 0 to N/2 contain all unique information; bins N/2+1 to N-1 are mirror images. Displaying only the first N/2+1 bins is called the "one-sided spectrum." It avoids confusing the user with duplicate information.

---

### Q8: How do you normalize the magnitude spectrum?

**A**: For a sinusoid `A × sin(2πfn/N)`, the FFT gives magnitude A×N/2 at bins +k and -k (both sides). To recover amplitude A from the one-sided spectrum:

- Non-DC, non-Nyquist bins: `normalized_magnitude = 2 × |X[k]| / N`
- DC (k=0) and Nyquist (k=N/2): `normalized_magnitude = |X[k]| / N`

The factor of 2 accounts for folding the two-sided spectrum into one side.

---

### Q9: What is spectral leakage and how can it be reduced?

**A**: Spectral leakage occurs when a signal's frequency falls between two FFT bins. Instead of a single sharp peak, energy "leaks" into neighboring bins, making frequency estimation inaccurate.

**Cause**: The FFT implicitly assumes the signal repeats periodically. If frequency isn't bin-aligned, the assumed periodicity causes discontinuities.

**Solution**: Apply a **window function** (e.g., Hann, Blackman-Harris) to the signal before FFT. Windows taper the signal to zero at edges, reducing discontinuities and leakage — at the cost of slightly lower frequency resolution.

---

### Q10: Why did you choose `std::complex<double>` instead of two separate arrays for real and imaginary parts?

**A**: `std::complex<double>` keeps the real and imaginary parts together, making the butterfly formula readable and matching the mathematical notation exactly. For maximum performance, a "split-complex" format (separate arrays for real and imaginary) can be better for SIMD vectorization. But for clarity and correctness in this project, `std::complex<double>` is the right choice.

---

### Q11: How would you adapt this for a fixed-point embedded system?

**A**: Replace `double` with `int16_t` (Q15 format, where 32767 ≈ 1.0). The butterfly multiplication becomes a fixed-point multiply-accumulate. The twiddle factors are stored as a pre-computed table in ROM (flash). This eliminates floating-point and avoids runtime `cos()`/`sin()` calls. The ARM CMSIS-DSP library (`arm_cfft_q15()`) provides a reference implementation.

---

### Q12: Why is FFT critical to OFDM in 4G/5G modems?

**A**: OFDM modulates data onto orthogonal subcarriers in the frequency domain. The transmitter uses IFFT to convert frequency-domain subcarrier symbols into a time-domain baseband waveform. The receiver uses FFT to reverse this and recover the subcarrier data. Without the FFT's O(N log N) efficiency, this would be computationally impossible in real time. A 4G LTE base station runs a 2048-point IFFT roughly 14,000 times per second.

---

### Q13: What is the space complexity of your FFT?

**A**: O(N) — in-place. The algorithm operates entirely on the input array, swapping elements and writing butterfly outputs back. No additional O(N) auxiliary array is needed (unlike a recursive out-of-place FFT). The only extra memory is O(log N) for loop variables, which is negligible.

---

### Q14: What would happen if you provided a non-power-of-2 input?

**A**: The bit-reversal permutation and butterfly pass structure are both based on `N = 2^k`. With a non-power-of-2 N, butterflies would either read/write out-of-bounds memory (undefined behavior) or produce incorrect results. My implementation validates N at the start of `fft()` and throws `std::invalid_argument` with a descriptive message.

---

### Q15: Why does the DFT of an impulse signal equal all ones?

**A**: The DFT formula is `X[k] = Σ x[n] × e^(-j·2πkn/N)`. For an impulse at n=0, x[0]=1 and all other x[n]=0. So X[k] = x[0] × e^0 = 1 × 1 = 1 for every k. The impulse has equal energy at all frequencies — it's a **flat spectrum**.

---

### Q16: How do you measure FFT execution time, and why is it important?

**A**: I use `std::chrono::high_resolution_clock`:

```cpp
auto t0 = std::chrono::high_resolution_clock::now();
fft::fft(data);
auto t1 = std::chrono::high_resolution_clock::now();
auto us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
```

In modem firmware, FFT latency directly impacts processing frame deadlines. LTE has strict 1 ms subframe timing. If FFT takes too long, the receiver misses its processing window.

---

### Q17: Your project uses `double` precision. In a real modem, would you use `float` or fixed-point?

**A**: Modern modem DSPs use **fixed-point 16-bit or 32-bit** arithmetic (Q15/Q31 format) for maximum throughput, low power, and small silicon area. Some high-performance modems use **single-precision float (`float`)** with SIMD. `double` is generally too slow and power-hungry for the critical FFT path in a mobile device. I used `double` in this project for numerical clarity and ease of verification.

---

### Q18: What is your peak detection algorithm and what are its weaknesses?

**A**: My algorithm finds local maxima in the magnitude spectrum — bins where the magnitude exceeds both neighbors and a minimum threshold (1% of the maximum magnitude). Then I sort by amplitude and return the top two.

**Weaknesses**:
1. Fails with spectral leakage (peak spreads across multiple bins).
2. Fails with noise (noise peaks can exceed signal peaks).
3. Cannot resolve two closely spaced frequencies within one bin's resolution.
4. Hardcoded 2-peak limit.

**Robust alternatives**: CFAR detection, matched filtering, MUSIC algorithm.

---

### Q19: Why is the Nyquist theorem important in your design?

**A**: The Nyquist theorem states that to faithfully digitize a signal of maximum frequency `f_max`, you need `Fs > 2 × f_max`. My project validates this: if the user passes frequencies above Fs/2, the program rejects them with an error.

**Consequence of violation**: If f2 = 600 Hz with Fs = 1024 Hz, the ADC would sample too slowly, and the 600 Hz signal would appear as a false frequency 1024 - 600 = 424 Hz (aliasing). This is why anti-aliasing filters are placed before every ADC in real systems.

---

### Q20: If an interviewer asked you to extend this to a real-time streaming FFT, how would you approach it?

**A**: I would implement **overlapping block processing (Short-Time Fourier Transform, STFT)**:

1. **Input buffer**: Ring buffer receiving ADC samples via interrupt/DMA.
2. **Block extraction**: When N new samples arrive, extract a block (with 50% overlap for STFT).
3. **Windowing**: Apply a Hann window to each block before FFT.
4. **FFT**: Run the N-point FFT on each block.
5. **Post-processing**: Compute magnitude, detect peaks, output to display or radio chain.
6. **Threading**: Use a producer-consumer model with a lockless FIFO — ADC fills the buffer, processing thread drains it.

For an embedded system: use DMA double-buffering (ping-pong buffers), so one buffer fills while the other is being processed. Triggered by a timer ISR at the ADC sample rate.

---

*End of INTERVIEW_PREP.md*
