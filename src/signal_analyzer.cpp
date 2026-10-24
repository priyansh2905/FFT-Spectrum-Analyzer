/**
 * @file signal_analyzer.cpp
 * @brief Signal generation and spectrum analysis — Implementation
 */

#include "signal_analyzer.h"
#include "fft.h"

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

static constexpr double PI = 3.14159265358979323846;

namespace analyzer {

// ============================================================
// validateParams
// ============================================================
void validateParams(const SignalParams& params) {
    // Check N is nonzero and power of two
    if (params.N == 0 || !fft::isPowerOfTwo(params.N)) {
        throw std::invalid_argument(
            "N must be a positive power of 2. Got: " + std::to_string(params.N)
        );
    }

    // Check frequencies are positive
    if (params.f1 <= 0.0 || params.f2 <= 0.0) {
        throw std::invalid_argument(
            "Signal frequencies must be positive. f1=" +
            std::to_string(params.f1) + " f2=" + std::to_string(params.f2)
        );
    }

    // Check sampling frequency is positive
    if (params.fs <= 0.0) {
        throw std::invalid_argument(
            "Sampling frequency must be positive. Fs=" + std::to_string(params.fs)
        );
    }

    // Nyquist criterion: Fs > 2 * f_max (strict inequality for aliasing safety)
    const double f_max = std::max(params.f1, params.f2);
    if (params.fs <= 2.0 * f_max) {
        throw std::invalid_argument(
            "Nyquist criterion violated: Fs (" + std::to_string(params.fs) +
            ") must be greater than 2 * f_max (" + std::to_string(2.0 * f_max) +
            "). Reduce signal frequencies or increase Fs."
        );
    }

    // Sanity: amplitudes should be positive
    if (params.amp1 <= 0.0 || params.amp2 <= 0.0) {
        throw std::invalid_argument("Amplitudes must be positive.");
    }
}

// ============================================================
// generateSignal
// ============================================================
ComplexVec generateSignal(const SignalParams& params) {
    // Allocate N complex samples; imaginary parts stay 0.0 for real signals
    ComplexVec signal(params.N);

    for (std::size_t n = 0; n < params.N; ++n) {
        // Time at sample n: t = n / Fs (in seconds)
        // x[n] = A1 * sin(2π * f1 * n/Fs) + A2 * sin(2π * f2 * n/Fs)
        const double t = static_cast<double>(n) / params.fs;

        const double sample =
            params.amp1 * std::sin(2.0 * PI * params.f1 * t) +
            params.amp2 * std::sin(2.0 * PI * params.f2 * t);

        signal[n] = Complex(sample, 0.0);  // Real-valued signal, imaginary = 0
    }

    return signal;
}

// ============================================================
// computeMagnitudeSpectrum
// ============================================================
std::vector<double> computeMagnitudeSpectrum(const ComplexVec& fftOutput) {
    const std::size_t N = fftOutput.size();
    // One-sided spectrum: bins 0 to N/2 inclusive → N/2 + 1 values
    const std::size_t numBins = N / 2 + 1;

    std::vector<double> magnitudes(numBins);

    for (std::size_t k = 0; k < numBins; ++k) {
        // Raw magnitude = |X[k]| = sqrt(Re²+ Im²)
        double mag = std::abs(fftOutput[k]);

        // Normalize to get true amplitude of sinusoid in the signal:
        //   - DC bin (k=0) and Nyquist bin (k=N/2): divide by N
        //   - All other bins: divide by N/2 (factor of 2 for one-sided spectrum)
        //
        // WHY? The FFT output X[k] for a sinusoid A*sin(2πfn/N) has
        // magnitude A*N/2 at ±k bins (both sides). For the one-sided
        // spectrum, we multiply by 2 and divide by N → amplitude = A.
        if (k == 0 || k == N / 2) {
            magnitudes[k] = mag / static_cast<double>(N);
        } else {
            magnitudes[k] = (2.0 * mag) / static_cast<double>(N);
        }
    }

    return magnitudes;
}

// ============================================================
// detectPeaks
// ============================================================
std::vector<FrequencyPeak> detectPeaks(
    const std::vector<double>& magnitudes,
    const SignalParams& params,
    std::size_t numPeaks
) {
    const std::size_t numBins = magnitudes.size();

    // Minimum threshold: ignore bins with amplitude below 1% of max
    const double maxMag = *std::max_element(magnitudes.begin(), magnitudes.end());
    const double threshold = 0.01 * maxMag;

    // Collect all local maxima: a bin is a peak if it is greater than
    // both its neighbours and above the threshold.
    std::vector<FrequencyPeak> peaks;

    for (std::size_t k = 1; k + 1 < numBins; ++k) {
        if (magnitudes[k] > magnitudes[k - 1] &&
            magnitudes[k] > magnitudes[k + 1] &&
            magnitudes[k] > threshold)
        {
            // Convert bin index to frequency: freq = k * Fs / N
            const double freqHz =
                static_cast<double>(k) * params.fs / static_cast<double>(params.N);

            peaks.push_back({freqHz, magnitudes[k], k});
        }
    }

    // Sort peaks by amplitude, descending
    std::sort(peaks.begin(), peaks.end(),
        [](const FrequencyPeak& a, const FrequencyPeak& b) {
            return a.normalizedAmp > b.normalizedAmp;
        }
    );

    // Return only the top numPeaks peaks
    if (peaks.size() > numPeaks) {
        peaks.resize(numPeaks);
    }

    // Sort by frequency for display (ascending)
    std::sort(peaks.begin(), peaks.end(),
        [](const FrequencyPeak& a, const FrequencyPeak& b) {
            return a.frequencyHz < b.frequencyHz;
        }
    );

    return peaks;
}

// ============================================================
// exportSignalToCSV
// ============================================================
void exportSignalToCSV(
    const ComplexVec& signal,
    const SignalParams& params,
    const std::string& filepath
) {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    file << "time_s,amplitude\n";
    file.precision(10);

    for (std::size_t n = 0; n < signal.size(); ++n) {
        const double t = static_cast<double>(n) / params.fs;
        file << t << "," << signal[n].real() << "\n";
    }

    file.close();
}

// ============================================================
// exportSpectrumToCSV
// ============================================================
void exportSpectrumToCSV(
    const std::vector<double>& magnitudes,
    const SignalParams& params,
    const std::string& filepath
) {
    std::ofstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    file << "frequency_hz,magnitude\n";
    file.precision(10);

    for (std::size_t k = 0; k < magnitudes.size(); ++k) {
        const double freqHz =
            static_cast<double>(k) * params.fs / static_cast<double>(params.N);
        file << freqHz << "," << magnitudes[k] << "\n";
    }

    file.close();
}

} // namespace analyzer
