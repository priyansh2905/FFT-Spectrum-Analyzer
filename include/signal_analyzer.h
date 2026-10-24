/**
 * @file signal_analyzer.h
 * @brief Signal generation and spectrum analysis — Declarations
 *
 * This module handles:
 *   1. Generating a test signal as a sum of sinusoids.
 *   2. Computing FFT magnitudes and converting bin indices to Hz.
 *   3. Detecting dominant frequency peaks in the spectrum.
 *   4. Exporting results to CSV for external visualization.
 *
 * SIGNAL MODEL:
 *   x[n] = A1 * sin(2π * f1 * n / Fs) + A2 * sin(2π * f2 * n / Fs)
 *
 * FREQUENCY BIN TO HZ CONVERSION:
 *   frequency_hz = k * Fs / N
 *   where k = bin index, Fs = sampling frequency, N = total samples
 *
 * AMPLITUDE NORMALIZATION (one-sided spectrum):
 *   For a real-valued signal, the two-sided FFT magnitude is
 *   symmetric. The one-sided spectrum magnitude is normalized by:
 *     - DC (k=0) and Nyquist (k=N/2): magnitude / N
 *     - All other bins:                 magnitude * 2 / N
 *
 * LIMITATION NOTE:
 *   Peak detection here finds the two largest magnitude peaks.
 *   This simple method works well for noiseless, bin-aligned
 *   sinusoids. With noise or spectral leakage (frequencies not
 *   exactly on a bin), this approach becomes unreliable.
 *   Production systems use windowing, interpolation, or MUSIC/ESPRIT
 *   algorithms for robust frequency estimation.
 *
 * NYQUIST THEOREM:
 *   To represent a signal of frequency f without aliasing, the
 *   sampling rate Fs must satisfy: Fs > 2 * f_max.
 *   For this project: Fs = 1024 Hz, f_max = 120 Hz → OK (512 > 120).
 */

#pragma once

#include "fft.h"
#include <vector>
#include <string>

namespace analyzer {

/**
 * @brief Parameters for signal generation and analysis.
 */
struct SignalParams {
    double fs;          ///< Sampling frequency in Hz
    std::size_t N;      ///< Number of samples (must be power of 2)
    double f1;          ///< First sinusoid frequency in Hz
    double f2;          ///< Second sinusoid frequency in Hz
    double amp1;        ///< Amplitude of first sinusoid (default 1.0)
    double amp2;        ///< Amplitude of second sinusoid (default 0.5)
};

/**
 * @brief Represents a detected frequency peak.
 */
struct FrequencyPeak {
    double frequencyHz;     ///< Detected frequency in Hz
    double normalizedAmp;   ///< Normalized amplitude (0.0 to 1.0)
    std::size_t binIndex;   ///< FFT bin index corresponding to this peak
};

/**
 * @brief Generates a time-domain signal as a sum of two sinusoids.
 *
 * x[n] = amp1 * sin(2π * f1 * n / Fs) + amp2 * sin(2π * f2 * n / Fs)
 *
 * @param params Signal parameters (Fs, N, f1, f2, amp1, amp2).
 * @return Vector of N complex values (imaginary parts = 0.0).
 */
ComplexVec generateSignal(const SignalParams& params);

/**
 * @brief Computes magnitude spectrum from FFT output.
 *
 * Returns |X[k]| for each bin k = 0 to N/2 (positive frequencies only).
 * Applies one-sided normalization as described above.
 *
 * @param fftOutput FFT result (N complex values).
 * @return Vector of N/2+1 normalized magnitude values.
 */
std::vector<double> computeMagnitudeSpectrum(const ComplexVec& fftOutput);

/**
 * @brief Detects the top 'numPeaks' dominant frequency peaks.
 *
 * Searches the magnitude spectrum for the largest distinct peaks.
 * Simple local-maximum approach: a bin is a peak if its magnitude
 * exceeds both neighbors and is above a minimum threshold.
 *
 * @param magnitudes  Normalized magnitude spectrum (N/2+1 values).
 * @param params      Signal parameters (used to convert bins to Hz).
 * @param numPeaks    Number of peaks to find (default = 2).
 * @return Sorted vector of FrequencyPeak (descending by amplitude).
 */
std::vector<FrequencyPeak> detectPeaks(
    const std::vector<double>& magnitudes,
    const SignalParams& params,
    std::size_t numPeaks = 2
);

/**
 * @brief Exports the time-domain signal to a CSV file.
 *
 * Format: two columns — time (seconds), amplitude (real part).
 *
 * @param signal  The time-domain signal samples.
 * @param params  Signal parameters (to compute time axis).
 * @param filepath Output file path.
 */
void exportSignalToCSV(
    const ComplexVec& signal,
    const SignalParams& params,
    const std::string& filepath
);

/**
 * @brief Exports the magnitude spectrum to a CSV file.
 *
 * Format: two columns — frequency (Hz), normalized amplitude.
 *
 * @param magnitudes  Normalized magnitude spectrum.
 * @param params      Signal parameters (to compute frequency axis).
 * @param filepath    Output file path.
 */
void exportSpectrumToCSV(
    const std::vector<double>& magnitudes,
    const SignalParams& params,
    const std::string& filepath
);

/**
 * @brief Validates signal parameters against Nyquist theorem.
 *
 * Checks:
 *   - N is a power of 2
 *   - Fs > 2 * max(f1, f2) (Nyquist criterion)
 *   - f1, f2 > 0
 *   - N > 0
 *
 * @param params Signal parameters to validate.
 * @throws std::invalid_argument if any parameter is invalid.
 */
void validateParams(const SignalParams& params);

} // namespace analyzer
