/**
 * @file main.cpp
 * @brief FFT-Based Digital Signal Spectrum Analyzer — Entry Point
 *
 * PROGRAM FLOW:
 * =============
 * 1. Parse command-line arguments (or use defaults).
 * 2. Validate signal parameters (Nyquist criterion, power-of-2, etc.).
 * 3. Generate time-domain signal: x[n] = A1*sin(2πf1*n/Fs) + A2*sin(2πf2*n/Fs).
 * 4. Run iterative Cooley-Tukey FFT (in-place, O(N log N)).
 * 5. Measure FFT execution time using std::chrono.
 * 6. Compute magnitude spectrum and normalize for one-sided display.
 * 7. Detect dominant frequency peaks.
 * 8. Display results to console.
 * 9. Export time-domain signal and spectrum to CSV files.
 * 10. (Optional) Compare FFT result with direct DFT on a small sub-array.
 *
 * USAGE:
 *   spectrum_analyzer [Fs] [N] [f1] [f2] [amp1] [amp2]
 *
 * DEFAULTS:
 *   Fs   = 1024 Hz    (Sampling frequency)
 *   N    = 1024       (Number of samples)
 *   f1   = 50 Hz      (First sinusoid frequency)
 *   f2   = 120 Hz     (Second sinusoid frequency)
 *   amp1 = 1.0        (Amplitude of first sinusoid)
 *   amp2 = 0.5        (Amplitude of second sinusoid)
 *
 * EXAMPLES:
 *   .\spectrum_analyzer.exe
 *   .\spectrum_analyzer.exe 2048 2048 100 300 1.0 0.8
 */

#include "fft.h"
#include "signal_analyzer.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <cmath>
#include <stdexcept>

// ============================================================
// printBanner
// ============================================================
static void printBanner() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════╗\n";
    std::cout << "║   FFT-Based Digital Signal Spectrum Analyzer  ║\n";
    std::cout << "║   Built with C++17 | Cooley-Tukey Radix-2    ║\n";
    std::cout << "╚══════════════════════════════════════════════╝\n";
    std::cout << "\n";
}

// ============================================================
// parseArgs — Parse command-line arguments into SignalParams
// ============================================================
static analyzer::SignalParams parseArgs(int argc, char* argv[]) {
    // Default parameters
    analyzer::SignalParams p;
    p.fs   = 1024.0;
    p.N    = 1024;
    p.f1   = 50.0;
    p.f2   = 120.0;
    p.amp1 = 1.0;
    p.amp2 = 0.5;

    // Override defaults if arguments are provided
    // Arguments: [Fs] [N] [f1] [f2] [amp1] [amp2]
    if (argc > 1) p.fs   = std::stod(argv[1]);
    if (argc > 2) p.N    = static_cast<std::size_t>(std::stoul(argv[2]));
    if (argc > 3) p.f1   = std::stod(argv[3]);
    if (argc > 4) p.f2   = std::stod(argv[4]);
    if (argc > 5) p.amp1 = std::stod(argv[5]);
    if (argc > 6) p.amp2 = std::stod(argv[6]);

    return p;
}

// ============================================================
// runDFTComparison — Compare FFT vs. DFT on a small subset
// ============================================================
static void runDFTComparison(const ComplexVec& originalSignal) {
    // Only run comparison on a small slice (first 64 samples)
    // to keep the O(N²) DFT fast.
    constexpr std::size_t COMPARE_N = 64;
    if (originalSignal.size() < COMPARE_N) return;

    std::cout << "\n[DFT vs FFT Comparison on first " << COMPARE_N << " samples]\n";
    std::cout << "  Computing DFT (O(N²))... ";

    // Take first 64 samples
    ComplexVec small(originalSignal.begin(), originalSignal.begin() + COMPARE_N);

    // Compute DFT directly
    auto t0 = std::chrono::high_resolution_clock::now();
    ComplexVec dftResult = fft::dft(small);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto dftUs = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

    // Compute FFT on the same data
    ComplexVec fftSmall = small;
    auto t2 = std::chrono::high_resolution_clock::now();
    fft::fft(fftSmall);
    auto t3 = std::chrono::high_resolution_clock::now();
    auto fftUs = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();

    // Compare: check max absolute difference across all bins
    double maxDiff = 0.0;
    for (std::size_t k = 0; k < COMPARE_N; ++k) {
        double diff = std::abs(dftResult[k] - fftSmall[k]);
        if (diff > maxDiff) maxDiff = diff;
    }

    std::cout << "Done.\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  DFT time  : " << dftUs << " µs  (O(N²), N=" << COMPARE_N << ")\n";
    std::cout << "  FFT time  : " << fftUs << " µs  (O(N log N), N=" << COMPARE_N << ")\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Max |DFT[k] - FFT[k]| = " << maxDiff
              << "  (should be near machine epsilon ~1e-13)\n";

    if (maxDiff < 1e-9) {
        std::cout << "  ✓ FFT output matches DFT output (numerical tolerance OK)\n";
    } else {
        std::cout << "  ✗ WARNING: FFT and DFT outputs differ more than expected!\n";
    }
    std::cout << "\n";
}

// ============================================================
// main
// ============================================================
int main(int argc, char* argv[]) {
    printBanner();

    // ── 1. Parse Parameters ──────────────────────────────────
    analyzer::SignalParams params;
    try {
        params = parseArgs(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << "\n";
        return 1;
    }

    // ── 2. Display Configuration ─────────────────────────────
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "Configuration\n";
    std::cout << "─────────────────────────────────────\n";
    std::cout << "  Sampling Frequency (Fs) : " << params.fs   << " Hz\n";
    std::cout << "  Number of Samples  (N)  : " << params.N    << "\n";
    std::cout << "  Signal Frequency 1 (f1) : " << params.f1   << " Hz  (amplitude " << params.amp1 << ")\n";
    std::cout << "  Signal Frequency 2 (f2) : " << params.f2   << " Hz  (amplitude " << params.amp2 << ")\n";
    std::cout << "  Frequency Resolution    : " << (params.fs / params.N) << " Hz/bin\n";
    std::cout << "  Nyquist Frequency       : " << (params.fs / 2.0) << " Hz\n";
    std::cout << "─────────────────────────────────────\n\n";

    // ── 3. Validate Parameters ───────────────────────────────
    try {
        analyzer::validateParams(params);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Parameter error: " << e.what() << "\n";
        return 1;
    }

    // ── 4. Generate Signal ───────────────────────────────────
    std::cout << "Generating signal...\n";
    const ComplexVec signal = analyzer::generateSignal(params);
    std::cout << "  Generated " << signal.size() << " samples.\n\n";

    // ── 5. Perform FFT ───────────────────────────────────────
    std::cout << "Performing FFT...\n";

    // Copy signal (FFT operates in-place; we need original for CSV export)
    ComplexVec fftData = signal;

    // Time the FFT execution
    const auto fftStart = std::chrono::high_resolution_clock::now();
    fft::fft(fftData);
    const auto fftEnd = std::chrono::high_resolution_clock::now();

    const auto fftDuration =
        std::chrono::duration_cast<std::chrono::microseconds>(fftEnd - fftStart);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  FFT execution time: " << fftDuration.count() << " µs\n";

    // Theoretical complexity note
    const std::size_t N = params.N;
    const double nLogN = static_cast<double>(N) * std::log2(static_cast<double>(N));
    const double nSquared = static_cast<double>(N) * static_cast<double>(N);
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "  Theoretical ops: FFT O(N log N) ≈ " << nLogN
              << " | DFT O(N²) ≈ " << nSquared << "\n\n";

    // ── 6. Compute Magnitude Spectrum ────────────────────────
    const std::vector<double> magnitudes = analyzer::computeMagnitudeSpectrum(fftData);

    // ── 7. Detect Peaks ──────────────────────────────────────
    const std::vector<analyzer::FrequencyPeak> peaks =
        analyzer::detectPeaks(magnitudes, params, 2);

    // ── 8. Display Results ───────────────────────────────────
    std::cout << "Detected Frequencies\n";
    std::cout << "─────────────────────────────────────\n";
    if (peaks.empty()) {
        std::cout << "  No significant peaks found.\n";
    } else {
        for (std::size_t i = 0; i < peaks.size(); ++i) {
            std::cout << std::fixed << std::setprecision(1);
            std::cout << "  Peak " << (i + 1)
                      << ": " << peaks[i].frequencyHz << " Hz"
                      << "  (normalized amplitude = ";
            std::cout << std::fixed << std::setprecision(4);
            std::cout << peaks[i].normalizedAmp
                      << ",  bin index = " << peaks[i].binIndex << ")\n";
        }
    }
    std::cout << "─────────────────────────────────────\n\n";

    // ── 9. DFT vs FFT Comparison ─────────────────────────────
    runDFTComparison(signal);

    // ── 10. Export CSV Files ─────────────────────────────────
    const std::string signalFile   = "results/signal.csv";
    const std::string spectrumFile = "results/spectrum.csv";

    std::cout << "Exporting results...\n";
    try {
        analyzer::exportSignalToCSV(signal, params, signalFile);
        std::cout << "  Signal   → " << signalFile << "\n";

        analyzer::exportSpectrumToCSV(magnitudes, params, spectrumFile);
        std::cout << "  Spectrum → " << spectrumFile << "\n";
    } catch (const std::runtime_error& e) {
        std::cerr << "  Warning: " << e.what() << "\n";
        std::cerr << "  (Create the 'results/' directory and re-run.)\n";
    }

    // ── 11. Summary ───────────────────────────────────────────
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════╗\n";
    std::cout << "║               FFT COMPLETE ✓                 ║\n";
    std::cout << "╚══════════════════════════════════════════════╝\n\n";

    std::cout << "To visualize results, run:\n";
    std::cout << "  python scripts\\plot.py\n\n";

    return 0;
}
