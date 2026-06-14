import numpy as np
import scipy.signal as sig
import pytest
import cmst_window as cmst


FILTER_POWER = 4
TAPS = 301
FS = 4096.0
N = 4096
BANDWIDTH = 200

# --- 1. THE ZERO-PHASE (SYMMETRY) TEST ---
def test_zero_phase_preservation():
    """An impulse must remain perfectly centered after filtering."""
    fs = FS
    n = N
    
    # Create an array of zeros with a single massive spike exactly in the middle
    impulse = np.zeros(n)
    center_idx = n // 2
    impulse[center_idx] = 1.0
    
    # Filter it
    filtered = cmst.cmst_lowpass(impulse, fs, f_high=BANDWIDTH, freq_power=2)
    
    # Find the new peak
    new_peak_idx = np.argmax(filtered)
    
    # Assert the peak has not shifted by a single index
    assert new_peak_idx == center_idx, "Filter injected a phase delay!"

# --- 2. THE AMPLITUDE CUTOFF TEST ---
def test_exact_50_percent_cutoff():
    """A sine wave exactly at f_high must be reduced to exactly 50% amplitude."""
    fs = FS
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
    fs = FS
    bandwidth = BANDWIDTH
    
    # Request an explicitly even number of taps
    coeffs = cmst.generate_cmst_fir(taps=TAPS+1, fs=fs, bandwidth=BANDWIDTH,freq_power=FILTER_POWER)
    
    # Assert the function caught it and appended the center tap
    assert len(coeffs) == TAPS+2, "FIR generator failed to force odd-symmetry!"


# --- 6. THE FREQUENCY RESPONSE SHAPE TEST ---
def test_fir_frequency_response_shape():
    """
    Calculates the exact transfer function of the FIR coefficients to verify
    the passband is preserved, the 50% cutoff is close, and the stopband is crushed.
    """
    fs = FS
    bandwidth = BANDWIDTH

    coeffs = cmst.generate_cmst_fir(taps=TAPS, fs=fs, bandwidth=bandwidth, freq_power=FILTER_POWER)

    # Calculate the frequency response
    w, h = sig.freqz(coeffs, worN=N, fs=fs)
    amplitude = np.abs(h)

    # Helper function to find the amplitude at a specific target frequency
    def get_amp_at_target(target_f):
        idx = np.argmin(np.abs(w - target_f))
        return amplitude[idx]

    amp_pass = get_amp_at_target(bandwidth / 2.0)
    amp_cut = get_amp_at_target(bandwidth)
    amp_stop = get_amp_at_target(bandwidth * 2.0)

    # Passband Verification
    assert amp_pass > 0.95, f"Passband attenuated too early! Amplitude at BANDWIDTH/2 Hz: {amp_pass:.3f}"

    # Cutoff Verification
    np.testing.assert_allclose(amp_cut, 0.5, rtol=0.01,
                               err_msg=f"Filter missed the 50% mark! Amplitude at BANDWIDTH Hz: {amp_cut:.3f}")

    # Stopband Verification
    # We assert it must be strictly less than 0.05 (-26dB).
    assert amp_stop < 0.05, f"Stopband is leaking noise! Amplitude at 200Hz: {amp_stop:.3f}"

# --- 7. DC GAIN NULLING TEST ---
def test_hp_blocks_dc():
    """
    A high-pass filter must have a DC gain of exactly 0.0.
    A low-pass filter should maintain a DC gain close to 1.0.
    """
    fs = FS
    bandwidth = 200.0
    taps = TAPS

    lp_coeffs = cmst.generate_cmst_lp_fir(TAPS, fs, bandwidth,freq_power=FILTER_POWER)
    hp_coeffs = cmst.generate_cmst_hp_fir(TAPS, fs, bandwidth,freq_power=FILTER_POWER)

    lp_dc_gain = np.sum(lp_coeffs)
    hp_dc_gain = np.sum(hp_coeffs)

    # Low-pass should sum to roughly 1.0 (depending on your internal normalization)
    # We use a loose tolerance here just to ensure it's passing DC.
    assert lp_dc_gain > 0.5, f"Low-pass filter is blocking DC! Gain: {lp_dc_gain:.3f}"

    # High-pass MUST sum to 0.0 to perfectly block DC
    np.testing.assert_allclose(hp_dc_gain, 0.0, atol=1e-10,
                               err_msg=f"High-pass filter failed to block DC! Gain: {hp_dc_gain:.8f}")


# --- 8. PERFECT COMPLEMENT (RECONSTRUCTION) TEST ---
def test_lp_hp_complementary_reconstruction():
    """
    Adding the LP and HP coefficients together should perfectly reconstruct
    an all-pass impulse (a Dirac delta function).
    """
    fs = FS
    bandwidth = 200.0
    taps = TAPS
    center_idx = taps // 2

    lp_coeffs = cmst.generate_cmst_lp_fir(TAPS, fs, bandwidth,freq_power=FILTER_POWER)
    hp_coeffs = cmst.generate_cmst_hp_fir(TAPS, fs, bandwidth,freq_power=FILTER_POWER)

    # Combine them
    combined_coeffs = lp_coeffs + hp_coeffs

    # The sum should be 1.0 at the center
    np.testing.assert_allclose(combined_coeffs[center_idx], 1.0, atol=1e-10,
                               err_msg="Combined center tap does not equal 1.0")

    # The sum should be 0.0 everywhere else (destructive interference)
    combined_coeffs[center_idx] = 0.0  # Zero out the center to test the tails
    tail_energy = np.sum(np.abs(combined_coeffs))

    np.testing.assert_allclose(tail_energy, 0.0, atol=1e-10,
                               err_msg="LP and HP tails did not perfectly cancel out!")


# --- 9. HIGH-PASS FREQUENCY RESPONSE TEST ---
def test_hp_frequency_response_shape():
    """
    Calculates the exact transfer function of the HP FIR coefficients to verify
    the stopband (low freqs) is crushed and the passband (high freqs) is preserved.
    """
    fs = FS
    bandwidth = BANDWIDTH  # Cutoff at 200 Hz
    taps = 301

    hp_coeffs = cmst.generate_cmst_hp_fir(TAPS, fs, bandwidth, freq_power=FILTER_POWER)

    w, h = sig.freqz(hp_coeffs, worN=N, fs=fs)
    amplitude = np.abs(h)

    def get_amp_at_target(target_f):
        idx = np.argmin(np.abs(w - target_f))
        return amplitude[idx]

    # For a HP filter, frequencies BELOW bandwidth should be blocked,
    # and frequencies ABOVE should pass.

    amp_blocked = get_amp_at_target(bandwidth/2.0)  # Deep in the stopband
    amp_cut = get_amp_at_target(bandwidth)  # Exactly at cutoff
    amp_passed = get_amp_at_target(bandwidth * 5)  # Deep in the passband

    # Stopband Verification (Low frequencies should be dead)
    assert amp_blocked < 0.05, f"High-pass is leaking low frequencies! Amplitude at BANDWIDTH/2 Hz: {amp_blocked:.3f}"

    # Cutoff Verification (Should still cross near 0.5)
    np.testing.assert_allclose(amp_cut, 0.5, rtol=0.01,
                               err_msg=f"High-pass missed the 50% cutoff! Amplitude at BANDWIDTH * 2 Hz: {amp_cut:.3f}")
    # Passband Verification (High frequencies should pass near 1.0)
    assert amp_passed > 0.95, f"High-pass attenuated the target signal! Amplitude at BANDWIDTH 0Hz: {amp_passed:.3f}"


def test_lp_hp_signal_recovery_in_noise():
    """
    Injects a composite signal (LP target + HP target + Gaussian noise) into both filters.
    Uses time-domain orthogonal demodulation to measure the recovered amplitude
    of the target signal and the rejection of the unwanted signal.
    """
    fs = FS
    bandwidth = BANDWIDTH  # Or whatever your standard test bandwidth is
    duration = 2.0
    t = np.arange(0, duration, 1.0 / fs)

    # 1. Define Frequencies
    f_lp = bandwidth / 1.5  # 50 Hz (Target for LP, Reject for HP)
    f_hp = bandwidth * 1.5  # 200 Hz (Reject for LP, Target for HP)

    # Generate a textbook Blackman window FIR filter with exactly 301 taps
    # Note: firwin uses a normalized cutoff frequency (Nyquist = 1.0)
    nyquist = fs / 2.0
    normalized_cutoff = bandwidth / nyquist

    # Low-pass Blackman
    blackman_lp_coeffs = sig.firwin(TAPS, normalized_cutoff, window='blackman')

    # High-pass Hamming
    hamming_hp_coeffs = sig.firwin(TAPS, normalized_cutoff, window='hamming', pass_zero=False)


    # 2. Generate pure signals (Amplitude = 1.0)
    s_lp = np.sin(2 * np.pi * f_lp * t)
    s_hp = np.sin(2 * np.pi * f_hp * t)

    # 3. Add Gaussian white noise
    noise_std = 0.2  # Moderate noise
    rng = np.random.default_rng(seed=42)  # Lock the random state
    noise = rng.normal(0, noise_std, len(t))
    
    # 4. Create the messy composite input signal
    x = s_lp + s_hp + noise

    # 5. Generate both filters
    lp_coeffs = cmst.generate_cmst_lp_fir(TAPS, fs, bandwidth, freq_power=6)
    # lp_coeffs = blackman_lp_coeffs
    hp_coeffs = cmst.generate_cmst_hp_fir(TAPS, fs, bandwidth, freq_power=6)
    # hp_coeffs = hamming_hp_coeffs

    # 6. Apply filters
    # mode='same' keeps the array length identical to 't', but creates edge transients
    y_lp = np.convolve(x, lp_coeffs, mode='same')
    y_hp = np.convolve(x, hp_coeffs, mode='same')

    # 7. Mask out the filter edge transients (first and last 0.25 seconds)
    # FIR filters need time to "fill up" their buffer.
    mid_mask = (t > 0.25) & (t < 1.75)

    y_lp_mid = y_lp[mid_mask]
    y_hp_mid = y_hp[mid_mask]
    s_lp_ref = s_lp[mid_mask]
    s_hp_ref = s_hp[mid_mask]

    # --- Helper: Lock-in Amplitude Measurement ---
    # Multiplying by the reference sine wave and taking the mean isolates that exact frequency.
    # Multiplying the result by 2 restores the peak amplitude measurement.
    def measure_amplitude(signal, reference):
        return 2.0 * np.abs(np.mean(signal * reference))

    # --- A. Low-Pass Filter Evaluation ---
    lp_recovered_target = measure_amplitude(y_lp_mid, s_lp_ref)
    lp_leaked_reject = measure_amplitude(y_lp_mid, s_hp_ref)

    lp_target_db = 20.0 * np.log10(max(lp_recovered_target, 1e-20))
    lp_leaked_db = 20.0 * np.log10(max(lp_leaked_reject, 1e-20))

    print(
        f"\nLP Output -> Target ({f_lp}Hz) Amp: {lp_recovered_target:.3f} ({lp_target_db:+.1f} dB) | Leaked ({f_hp}Hz) Amp: {lp_leaked_reject:.4f} ({lp_leaked_db:.1f} dB)")


    # --- B. High-Pass Filter Evaluation ---
    hp_recovered_target = measure_amplitude(y_hp_mid, s_hp_ref)
    hp_leaked_reject = measure_amplitude(y_hp_mid, s_lp_ref)

    hp_target_db = 20.0 * np.log10(max(hp_recovered_target, 1e-100))
    hp_leaked_db = 20.0 * np.log10(max(hp_leaked_reject, 1e-100))

    print(
        f"HP Output -> Target ({f_hp}Hz) Amp: {hp_recovered_target:.3f} ({hp_target_db:+.1f} dB) | Leaked ({f_lp}Hz) Amp: {hp_leaked_reject:.4f} ({hp_leaked_db:.1f} dB)")


    assert lp_recovered_target > 0.95, f"LP failed to recover signal! Amp: {lp_recovered_target:.3f}"
    assert lp_leaked_reject < 0.05, f"LP failed to reject high frequencies! Leakage: {lp_leaked_reject:.3f}"

    assert hp_recovered_target > 0.95, f"HP failed to recover signal! Amp: {hp_recovered_target:.3f}"
    assert hp_leaked_reject < 0.05, f"HP failed to reject low frequencies! Leakage: {hp_leaked_reject:.3f}"