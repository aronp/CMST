import numpy as np
import scipy.signal as sig
import pytest
import cmst_window as cmst

# --- 1. THE ZERO-PHASE (SYMMETRY) TEST ---
def test_zero_phase_preservation():
    """An impulse must remain perfectly centered after filtering."""
    fs = 4096.0
    N = 4096
    
    # Create an array of zeros with a single massive spike exactly in the middle
    impulse = np.zeros(N)
    center_idx = N // 2
    impulse[center_idx] = 1.0
    
    # Filter it
    filtered = cmst.cmst_lowpass(impulse, fs, f_high=100.0, freq_power=2)
    
    # Find the new peak
    new_peak_idx = np.argmax(filtered)
    
    # Assert the peak has not shifted by a single index
    assert new_peak_idx == center_idx, "Filter injected a phase delay!"

# --- 2. THE AMPLITUDE CUTOFF TEST ---
def test_exact_50_percent_cutoff():
    """A sine wave exactly at f_high must be reduced to exactly 50% amplitude."""
    fs = 4096.0
    duration = 2.0
    t = np.arange(0, duration, 1.0/fs)
    f_target = 150.0
    
    # Create a pure sine wave at exactly the cutoff frequency
    sine_wave = np.sin(2 * np.pi * f_target * t)
    
    # Filter it using f_high = f_target
    filtered = cmst.cmst_lowpass(sine_wave, fs, f_high=f_target, freq_power=2)
    
    # Measure the RMS amplitude of the middle section (ignoring boundaries)
    mid_mask = (t > 0.5) & (t < 1.5)
    original_rms = np.sqrt(np.mean(sine_wave[mid_mask]**2))
    filtered_rms = np.sqrt(np.mean(filtered[mid_mask]**2))
    
    ratio = filtered_rms / original_rms
    
    # Assert the ratio is 0.5 (allowing for tiny floating-point rounding margins)
    np.testing.assert_allclose(ratio, 0.5, rtol=1e-3, 
                               err_msg="Cutoff amplitude is not exactly 50%")

def test_fir_generator_forces_odd_taps():
    """Requesting an even number of taps must return an odd-length array."""
    fs = 4096.0
    bandwidth = 100.0
    
    # Request an explicitly even number of taps (100)
    coeffs = cmst.generate_cmst_fir(taps=100, fs=fs, bandwidth=bandwidth)
    
    # Assert the function caught it and appended the center tap
    assert len(coeffs) == 101, "FIR generator failed to force odd-symmetry!"


# --- 6. THE FREQUENCY RESPONSE SHAPE TEST ---
def test_fir_frequency_response_shape():
    """
    Calculates the exact transfer function of the FIR coefficients to verify
    the passband is preserved, the 50% cutoff is close, and the stopband is crushed.
    """
    fs = 4096.0
    bandwidth = 100.0

    # Generate the 101-tap array
    coeffs = cmst.generate_cmst_fir(taps=101, fs=fs, bandwidth=bandwidth, freq_power=2)

    # Calculate the frequency response
    # worN=8192 provides a dense, high-resolution frequency axis to interrogate
    w, h = sig.freqz(coeffs, worN=8192, fs=fs)
    amplitude = np.abs(h)

    # Helper function to find the amplitude at a specific target frequency
    def get_amp_at_target(target_f):
        idx = np.argmin(np.abs(w - target_f))
        return amplitude[idx]

    amp_pass = get_amp_at_target(bandwidth / 2.0)  # 50 Hz
    amp_cut = get_amp_at_target(bandwidth)  # 100 Hz
    amp_stop = get_amp_at_target(bandwidth * 2.0)  # 200 Hz

    # Passband Verification
    assert amp_pass > 0.95, f"Passband attenuated too early! Amplitude at 50Hz: {amp_pass:.3f}"

    # Cutoff Verification
    np.testing.assert_allclose(amp_cut, 0.5, rtol=0.20,
                               err_msg=f"Filter missed the 50% mark! Amplitude at 100Hz: {amp_cut:.3f}")

    # Stopband Verification
    # The ripples for a 101-tap p=2 filter hover around -40dB (0.01 amplitude).
    # We assert it must be strictly less than 0.05 (-26dB).
    assert amp_stop < 0.05, f"Stopband is leaking noise! Amplitude at 200Hz: {amp_stop:.3f}"


