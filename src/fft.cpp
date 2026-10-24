/**
 * @file fft.cpp
 * @brief Iterative Cooley-Tukey Radix-2 FFT — Implementation
 *
 * ALGORITHM OVERVIEW:
 * ==================
 * The Cooley-Tukey FFT works by recursively splitting the N-point DFT
 * into smaller DFTs. The iterative version (bottom-up) avoids recursion
 * overhead and is memory-efficient — ideal for embedded systems.
 *
 * STEP 1: Bit-Reversal Permutation
 *   Before any computation, we rearrange the input array so that
 *   element at index i moves to index bit_reverse(i, log2(N)).
 *
 *   WHY? The recursive Cooley-Tukey algorithm splits the array into
 *   even and odd indices repeatedly. When unrolled, this produces
 *   a bit-reversed ordering. By pre-permuting the input, we can
 *   compute FFT in-place with a simple loop structure.
 *
 *   Example (N=8, 3 bits):
 *     Index 0 = 000 → 000 = 0  (stays)
 *     Index 1 = 001 → 100 = 4  (swapped)
 *     Index 2 = 010 → 010 = 2  (stays)
 *     Index 3 = 011 → 110 = 6  (swapped)
 *     Index 4 = 100 → 001 = 1  (swapped)
 *     Index 5 = 101 → 101 = 5  (stays)
 *     Index 6 = 110 → 011 = 3  (swapped)
 *     Index 7 = 111 → 111 = 7  (stays)
 *
 * STEP 2: Butterfly Passes (log₂N passes total)
 *   In each pass s = 1, 2, ..., log₂N:
 *     - We process groups of 2^s elements.
 *     - Each group gets log₂(2^s) = s butterfly operations.
 *     - Twiddle factor: W_m = e^(-j*2π/m), m = 2^s
 *
 *   Butterfly operation on element pair (u, t):
 *     t  = W * x[j + m/2]
 *     x[j + m/2] = x[j] - t   (upper butterfly output)
 *     x[j]       = x[j] + t   (lower butterfly output)
 *
 *   After log₂N passes, x contains the FFT output X[k].
 *
 * COMPLEXITY ANALYSIS:
 *   - log₂N butterfly passes
 *   - N/2 butterfly operations per pass
 *   - Each butterfly = O(1) (1 complex multiply + 2 adds)
 *   - Total: O(N/2 * log₂N) = O(N log N)
 *
 *   vs DFT: O(N²)
 *   For N=1024: FFT ≈ 5120 operations, DFT ≈ 1,048,576 operations
 *
 * SPACE COMPLEXITY: O(N) — in-place, no extra arrays needed.
 */

#include "fft.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

// π constant — more readable than M_PI (not guaranteed in C++17)
static constexpr double PI = 3.14159265358979323846;

namespace fft {

// ============================================================
// isPowerOfTwo
// ============================================================
bool isPowerOfTwo(std::size_t n) {
    // A power of two has exactly one bit set.
    // n & (n-1) clears the lowest set bit. If result is 0, n is a power of 2.
    // Edge case: n=0 is NOT a power of two.
    return (n > 0) && ((n & (n - 1)) == 0);
}

// ============================================================
// bitReversalPermutation
// ============================================================
void bitReversalPermutation(ComplexVec& x) {
    const std::size_t N = x.size();

    // Compute the number of bits needed to represent indices 0..N-1.
    // For N=1024, bits = 10.
    std::size_t bits = 0;
    std::size_t temp = N;
    while (temp > 1) {
        temp >>= 1;   // Right shift = divide by 2
        ++bits;
    }

    // For each index i, compute its bit-reversed index j.
    // Only swap when j > i to avoid double-swapping.
    for (std::size_t i = 0; i < N; ++i) {
        std::size_t j = 0;       // Will hold bit-reversed index
        std::size_t val = i;     // Current index to reverse

        for (std::size_t b = 0; b < bits; ++b) {
            // Take the least significant bit of val and append it to j
            j = (j << 1) | (val & 1);
            val >>= 1;
        }

        // Swap only once (when j > i, to avoid re-swapping)
        if (j > i) {
            std::swap(x[i], x[j]);
        }
    }
}

// ============================================================
// fft — Iterative Cooley-Tukey Radix-2 FFT
// ============================================================
void fft(ComplexVec& x) {
    const std::size_t N = x.size();

    // Validate: N must be a power of two
    if (!isPowerOfTwo(N)) {
        throw std::invalid_argument(
            "FFT input size must be a power of 2. Got: " + std::to_string(N)
        );
    }

    // Edge case: 1-point FFT is a no-op (X[0] = x[0])
    if (N <= 1) return;

    // Step 1: Bit-reversal permutation
    bitReversalPermutation(x);

    // Step 2: Butterfly passes
    // s = current stage (1 to log2(N))
    // m = size of sub-DFT being computed = 2^s
    for (std::size_t s = 1; (1u << s) <= N; ++s) {
        const std::size_t m = 1u << s;        // m = 2^s (sub-DFT size)
        const std::size_t half_m = m >> 1;    // m/2 (number of butterfly pairs per group)

        // Twiddle factor base for this stage: W_m = e^(-j*2π/m)
        // This is the principal m-th root of unity.
        const Complex W_m = Complex(
            std::cos(-2.0 * PI / static_cast<double>(m)),
            std::sin(-2.0 * PI / static_cast<double>(m))
        );

        // Process each group of m elements
        for (std::size_t k = 0; k < N; k += m) {
            Complex W = Complex(1.0, 0.0);  // Current twiddle: W_m^j, starts at W_m^0 = 1

            // Butterfly operations within this group
            for (std::size_t j = 0; j < half_m; ++j) {
                // Indices of the two elements in this butterfly
                const std::size_t idx_top = k + j;
                const std::size_t idx_bot = k + j + half_m;

                // Butterfly computation:
                //   t   = W * x[bottom]
                //   x[top]    = x[top] + t
                //   x[bottom] = x[top] - t
                const Complex t = W * x[idx_bot];
                x[idx_bot] = x[idx_top] - t;
                x[idx_top] = x[idx_top] + t;

                // Advance twiddle factor: W_m^(j+1) = W_m^j * W_m
                W *= W_m;
            }
        }
    }
    // After all stages, x now contains FFT output X[0..N-1]
}

// ============================================================
// ifft — Inverse FFT using conjugation method
// ============================================================
void ifft(ComplexVec& x) {
    const std::size_t N = x.size();

    if (!isPowerOfTwo(N)) {
        throw std::invalid_argument(
            "IFFT input size must be a power of 2. Got: " + std::to_string(N)
        );
    }

    // IFFT via conjugation trick:
    //   IFFT(X) = conj(FFT(conj(X))) / N
    //
    // WHY THIS WORKS:
    //   The IFFT formula is: x[n] = (1/N) * Σ X[k] * e^(+j*2πkn/N)
    //   The FFT formula is:  X[k] = Σ x[n] * e^(-j*2πkn/N)
    //   Taking conjugate flips the sign of the exponent, converting
    //   FFT into IFFT (except for the 1/N normalization factor).

    // Step 1: Conjugate all elements
    for (auto& val : x) {
        val = std::conj(val);
    }

    // Step 2: Apply FFT (forward FFT on conjugated input)
    fft(x);

    // Step 3: Conjugate again and divide by N
    const double inv_N = 1.0 / static_cast<double>(N);
    for (auto& val : x) {
        val = std::conj(val) * inv_N;
    }
}

// ============================================================
// dft — Direct O(N²) DFT for reference/verification
// ============================================================
ComplexVec dft(const ComplexVec& x) {
    const std::size_t N = x.size();
    ComplexVec X(N, Complex(0.0, 0.0));

    // Naive DFT: X[k] = Σ_{n=0}^{N-1} x[n] * e^(-j*2π*k*n/N)
    for (std::size_t k = 0; k < N; ++k) {
        for (std::size_t n = 0; n < N; ++n) {
            // Compute twiddle factor e^(-j*2π*k*n/N)
            const double angle = -2.0 * PI * static_cast<double>(k) *
                                  static_cast<double>(n) / static_cast<double>(N);
            X[k] += x[n] * Complex(std::cos(angle), std::sin(angle));
        }
    }

    return X;
}

} // namespace fft
