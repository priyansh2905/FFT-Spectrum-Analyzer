/**
 * @file fft.h
 * @brief FFT (Fast Fourier Transform) — Declarations
 *
 * ============================================================
 * BACKGROUND: What is DFT and FFT?
 * ============================================================
 *
 * DFT (Discrete Fourier Transform):
 *   Converts a time-domain signal of N samples into its
 *   frequency-domain representation (N complex values).
 *   Formula: X[k] = Σ x[n] * e^(-j*2π*k*n/N),  n = 0..N-1
 *
 *   Complexity: O(N²) — for each of N output bins, we sum N terms.
 *   For N = 1024: ~1,000,000 operations. Too slow for real-time.
 *
 * FFT (Fast Fourier Transform):
 *   An algorithm that computes the SAME result as DFT but in
 *   O(N log₂N) operations by exploiting symmetry and periodicity
 *   of the complex exponential (twiddle factors).
 *
 *   For N = 1024: ~10,240 operations. ~100× faster than DFT.
 *
 * This file implements the Cooley-Tukey Radix-2 FFT algorithm:
 *   The most common FFT algorithm, requires N = power of 2.
 *   It recursively splits the DFT of size N into two DFTs of
 *   size N/2, until trivial 1-point DFTs are reached.
 *
 * BUTTERFLY OPERATION:
 *   The core computation unit of FFT. Takes two complex numbers
 *   (a, b) and a twiddle factor W = e^(-j*2π*k/N):
 *       a' = a + W * b
 *       b' = a - W * b
 *   Called "butterfly" because the signal flow diagram looks
 *   like butterfly wings.
 *
 * BIT-REVERSAL PERMUTATION:
 *   The iterative FFT processes data in a specific order. Before
 *   the butterfly passes, input samples must be reordered by
 *   reversing the bits of each sample's index. For N=8, index 3
 *   (binary 011) goes to index 6 (binary 110).
 *   This ensures data is in the correct position for in-place
 *   butterfly computations to work correctly.
 *
 * RELEVANCE TO WIRELESS / OFDM:
 *   OFDM (Orthogonal Frequency Division Multiplexing), used in
 *   4G LTE, 5G NR, WiFi, and DVB, maps data onto orthogonal
 *   subcarriers using IFFT at the transmitter and FFT at the
 *   receiver. A 4G LTE modem runs IFFT/FFT thousands of times
 *   per second (e.g., 2048-point FFT at 30.72 MHz). This is why
 *   efficient FFT implementation is critical in cellular modem
 *   firmware and PHY layer hardware/software.
 * ============================================================
 */

#pragma once

#include <vector>
#include <complex>
#include <cstddef>

// Type alias for clarity
using Complex = std::complex<double>;
using ComplexVec = std::vector<Complex>;

namespace fft {

/**
 * @brief Checks if n is a power of two (n > 0).
 */
bool isPowerOfTwo(std::size_t n);

/**
 * @brief Computes the bit-reversal permutation of the input.
 *
 * Reorders x so that x[i] is placed at x[bit_reverse(i)].
 * This step is required before the iterative butterfly passes.
 *
 * @param x Input/output complex vector (modified in-place).
 */
void bitReversalPermutation(ComplexVec& x);

/**
 * @brief Performs in-place iterative Cooley-Tukey Radix-2 FFT.
 *
 * Computes the forward DFT: X[k] = Σ x[n] * e^(-j*2π*k*n/N)
 * Convention: unnormalized forward DFT (standard convention).
 *
 * @param x Complex vector of N samples (N must be power of 2).
 * @throws std::invalid_argument if N is not a power of two.
 */
void fft(ComplexVec& x);

/**
 * @brief Performs in-place Inverse FFT (IFFT).
 *
 * Uses the conjugation method:
 *   IFFT(x) = conj(FFT(conj(x))) / N
 *
 * This reconstructs the original time-domain signal from
 * frequency-domain data.
 *
 * @param x Complex frequency-domain vector (modified in-place).
 * @throws std::invalid_argument if N is not a power of two.
 */
void ifft(ComplexVec& x);

/**
 * @brief Computes the DFT directly (O(N²)) for reference/testing.
 *
 * Used to verify FFT correctness on small inputs (N ≤ 256).
 * Implements the naive O(N²) DFT formula without any optimization.
 *
 * @param x Input time-domain complex vector.
 * @return Complex vector of DFT output.
 */
ComplexVec dft(const ComplexVec& x);

} // namespace fft
