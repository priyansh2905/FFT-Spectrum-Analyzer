/**
 * @file test_fft.cpp
 * @brief Automated tests for FFT implementation.
 *
 * Tests:
 *   1. FFT of an impulse signal (δ[n]): output should be all 1.0+0j
 *   2. FFT of a constant signal (DC): output should be N at bin 0, 0 elsewhere
 *   3. FFT of a 50 Hz sinusoid: peak at bin 50 (for N=Fs=1024)
 *   4. IFFT(FFT(x)) = x (reconstruction within floating-point tolerance)
 *   5. Invalid input (N not power of 2): should throw std::invalid_argument
 *
 * Floating-point tolerance:
 *   Due to finite-precision arithmetic, computed values may differ
 *   from exact values by small amounts (~1e-12). We use a tolerance
 *   of 1e-9 to check approximate equality.
 *
 * COMPILATION (via CTest):
 *   cmake --build build --target test_fft
 *   ctest --test-dir build --verbose
 */

#include "fft.h"
#include "signal_analyzer.h"

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <string>
#include <cassert>

// ============================================================
// Test infrastructure
// ============================================================

static int tests_run    = 0;
static int tests_passed = 0;
static int tests_failed = 0;

// Tolerance for floating-point comparisons
static constexpr double TOL = 1e-9;

/**
 * @brief Check a condition and report pass/fail.
 */
static void check(bool condition, const std::string& testName) {
    ++tests_run;
    if (condition) {
        ++tests_passed;
        std::cout << "  [PASS] " << testName << "\n";
    } else {
        ++tests_failed;
        std::cout << "  [FAIL] " << testName << "\n";
    }
}

/**
 * @brief Check that |a - b| < tolerance.
 */
static bool approxEqual(double a, double b, double tol = TOL) {
    return std::abs(a - b) < tol;
}

/**
 * @brief Check that two complex numbers are approximately equal.
 */
static bool approxEqualC(const Complex& a, const Complex& b, double tol = TOL) {
    return std::abs(a - b) < tol;
}

// ============================================================
// Test 1: FFT of an Impulse Signal
// ============================================================
/**
 * Mathematical property:
 *   x[0] = 1, x[1..N-1] = 0  →  X[k] = 1  for all k
 *
 * The DFT of a unit impulse at n=0 is a flat spectrum.
 * This is a fundamental property of the Fourier transform.
 */
static void testImpulseFFT() {
    std::cout << "\n[Test 1: FFT of Impulse Signal]\n";
    constexpr std::size_t N = 8;

    ComplexVec x(N, Complex(0.0, 0.0));
    x[0] = Complex(1.0, 0.0);  // Impulse at n=0

    fft::fft(x);

    // All output bins should equal 1.0 + 0j
    bool allOnes = true;
    for (std::size_t k = 0; k < N; ++k) {
        if (!approxEqualC(x[k], Complex(1.0, 0.0))) {
            allOnes = false;
            std::cout << "  Bin " << k << ": " << x[k] << " (expected 1+0j)\n";
        }
    }
    check(allOnes, "Impulse FFT: all bins = 1+0j");
}

// ============================================================
// Test 2: FFT of a Constant (DC) Signal
// ============================================================
/**
 * Mathematical property:
 *   x[n] = C for all n  →  X[0] = N*C,  X[k] = 0 for k ≠ 0
 *
 * A constant signal has energy only at DC (0 Hz).
 */
static void testDCSignalFFT() {
    std::cout << "\n[Test 2: FFT of Constant (DC) Signal]\n";
    constexpr std::size_t N = 16;
    constexpr double C = 3.0;

    ComplexVec x(N, Complex(C, 0.0));  // All samples = C

    fft::fft(x);

    // X[0] should equal N * C
    check(approxEqualC(x[0], Complex(static_cast<double>(N) * C, 0.0)),
          "DC FFT: X[0] = N*C = " + std::to_string(N * C));

    // All other bins should be 0
    bool allZero = true;
    for (std::size_t k = 1; k < N; ++k) {
        if (!approxEqualC(x[k], Complex(0.0, 0.0))) {
            allZero = false;
            std::cout << "  Bin " << k << ": " << x[k] << " (expected 0+0j)\n";
        }
    }
    check(allZero, "DC FFT: X[k] = 0 for k != 0");
}

// ============================================================
// Test 3: FFT of a 50 Hz Sinusoid
// ============================================================
/**
 * With N = Fs = 1024, a 50 Hz sinusoid is bin-aligned.
 * The FFT should show a dominant peak at bin 50.
 */
static void testSinusoidDetection() {
    std::cout << "\n[Test 3: FFT of a 50 Hz Sinusoid (Fs=1024, N=1024)]\n";

    analyzer::SignalParams params;
    params.fs   = 1024.0;
    params.N    = 1024;
    params.f1   = 50.0;
    params.f2   = 50.0;   // Both frequencies the same (single sinusoid)
    params.amp1 = 1.0;
    params.amp2 = 0.0;    // Zero amplitude for second sinusoid

    // Generate x[n] = sin(2π * 50 * n / 1024)
    // Since amp2 = 0.0, this would fail validateParams for amp2 <= 0.
    // So we generate the signal manually here.
    ComplexVec signal(params.N);
    constexpr double PI = 3.14159265358979323846;
    for (std::size_t n = 0; n < params.N; ++n) {
        const double t = static_cast<double>(n) / params.fs;
        signal[n] = Complex(std::sin(2.0 * PI * 50.0 * t), 0.0);
    }

    fft::fft(signal);

    // Compute magnitudes for one-sided spectrum
    const std::size_t numBins = params.N / 2 + 1;
    std::size_t peakBin = 0;
    double peakMag = 0.0;

    for (std::size_t k = 0; k < numBins; ++k) {
        double mag = std::abs(signal[k]);
        if (mag > peakMag) {
            peakMag = mag;
            peakBin = k;
        }
    }

    // The dominant peak should be at bin 50 (= 50 Hz)
    check(peakBin == 50,
          "Sinusoid detection: peak at bin 50 (got bin " + std::to_string(peakBin) + ")");

    // Magnitude at bin 50 should ≈ N/2 = 512 (for A=1.0, unnormalized)
    check(approxEqual(peakMag, static_cast<double>(params.N) / 2.0, 1.0),
          "Sinusoid magnitude: |X[50]| ≈ N/2 = 512 (got " +
          std::to_string(static_cast<long long>(peakMag)) + ")");
}

// ============================================================
// Test 4: IFFT(FFT(x)) = x (Round-Trip Reconstruction)
// ============================================================
/**
 * Applying FFT then IFFT should recover the original signal
 * within floating-point rounding error.
 */
static void testIFFTReconstruction() {
    std::cout << "\n[Test 4: IFFT(FFT(x)) Round-Trip Reconstruction]\n";
    constexpr std::size_t N = 64;

    // Create a non-trivial test signal
    ComplexVec original(N);
    constexpr double PI = 3.14159265358979323846;
    for (std::size_t n = 0; n < N; ++n) {
        original[n] = Complex(
            std::sin(2.0 * PI * 5.0 * n / N) + 0.5 * std::cos(2.0 * PI * 13.0 * n / N),
            0.0
        );
    }

    // Apply FFT then IFFT
    ComplexVec x = original;
    fft::fft(x);
    fft::ifft(x);

    // Check that each sample matches the original within tolerance
    double maxError = 0.0;
    for (std::size_t n = 0; n < N; ++n) {
        double err = std::abs(x[n] - original[n]);
        if (err > maxError) maxError = err;
    }

    check(maxError < TOL,
          "IFFT reconstruction: max error = " + std::to_string(maxError) +
          " (threshold " + std::to_string(TOL) + ")");
}

// ============================================================
// Test 5: Invalid Input (Non-Power-of-Two N)
// ============================================================
/**
 * FFT should throw std::invalid_argument for N = 100 (not a power of 2).
 */
static void testInvalidInput() {
    std::cout << "\n[Test 5: Invalid Input — Non-Power-of-Two N]\n";

    ComplexVec x(100, Complex(1.0, 0.0));  // N=100, not a power of 2

    bool threw = false;
    try {
        fft::fft(x);
    } catch (const std::invalid_argument& e) {
        threw = true;
        std::cout << "  Caught expected exception: " << e.what() << "\n";
    } catch (...) {
        std::cout << "  Unexpected exception type caught.\n";
    }

    check(threw, "FFT throws invalid_argument for N=100");
}

// ============================================================
// Test 6: isPowerOfTwo utility
// ============================================================
static void testIsPowerOfTwo() {
    std::cout << "\n[Test 6: isPowerOfTwo utility]\n";
    check( fft::isPowerOfTwo(1),    "isPowerOfTwo(1) = true");
    check( fft::isPowerOfTwo(2),    "isPowerOfTwo(2) = true");
    check( fft::isPowerOfTwo(4),    "isPowerOfTwo(4) = true");
    check( fft::isPowerOfTwo(1024), "isPowerOfTwo(1024) = true");
    check(!fft::isPowerOfTwo(0),    "isPowerOfTwo(0) = false");
    check(!fft::isPowerOfTwo(3),    "isPowerOfTwo(3) = false");
    check(!fft::isPowerOfTwo(100),  "isPowerOfTwo(100) = false");
    check(!fft::isPowerOfTwo(1023), "isPowerOfTwo(1023) = false");
}

// ============================================================
// Test 7: Two-frequency signal peak detection
// ============================================================
/**
 * Full pipeline test: generate 50+120 Hz signal, run FFT,
 * detect peaks, verify both frequencies are found.
 */
static void testTwoFrequencyDetection() {
    std::cout << "\n[Test 7: Full Pipeline — 50 Hz + 120 Hz Detection]\n";

    analyzer::SignalParams params;
    params.fs   = 1024.0;
    params.N    = 1024;
    params.f1   = 50.0;
    params.f2   = 120.0;
    params.amp1 = 1.0;
    params.amp2 = 0.5;

    // Generate and process signal
    ComplexVec signal = analyzer::generateSignal(params);
    fft::fft(signal);
    std::vector<double> magnitudes = analyzer::computeMagnitudeSpectrum(signal);
    std::vector<analyzer::FrequencyPeak> peaks = analyzer::detectPeaks(magnitudes, params, 2);

    check(peaks.size() == 2, "Two peaks detected");

    if (peaks.size() == 2) {
        // Peaks are sorted by frequency ascending
        check(approxEqual(peaks[0].frequencyHz, 50.0, 0.5),
              "Peak 1 = 50 Hz (got " + std::to_string(peaks[0].frequencyHz) + " Hz)");
        check(approxEqual(peaks[1].frequencyHz, 120.0, 0.5),
              "Peak 2 = 120 Hz (got " + std::to_string(peaks[1].frequencyHz) + " Hz)");
    }
}

// ============================================================
// main
// ============================================================
int main() {
    std::cout << "======================================\n";
    std::cout << "  FFT Spectrum Analyzer — Test Suite  \n";
    std::cout << "======================================\n";

    testIsPowerOfTwo();
    testImpulseFFT();
    testDCSignalFFT();
    testSinusoidDetection();
    testIFFTReconstruction();
    testInvalidInput();
    testTwoFrequencyDetection();

    std::cout << "\n======================================\n";
    std::cout << "Results: " << tests_passed << "/" << tests_run
              << " tests passed";
    if (tests_failed > 0) {
        std::cout << "  (" << tests_failed << " FAILED)";
    }
    std::cout << "\n======================================\n";

    // Return 0 if all tests pass, 1 otherwise (CTest checks exit code)
    return (tests_failed == 0) ? 0 : 1;
}
