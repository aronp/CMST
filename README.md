# cmst-window
**The analytically sound, zero-preserving, interlace-preserving compact window.**
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Stability: Verified](https://img.shields.io/badge/Stability-Verified-green)](https://github.com/aronp/CMST)
[![Theory: CMST](https://img.shields.io/badge/Theory-CMST-purple)](https://github.com/yourusername/cmst-window)

### 🚀 The Problem
Signal processing engineers often need a "Flat-Top" window to preserve signal amplitude while filtering noise. Existing standard solutions force a dangerous compromise:

* **Truncated Super-Gaussians:** Create **geometric singularities** (infinite acceleration) at the boundaries, causing ringing in control loops.
* **Planck-taper (LIGO Standard):** Theoretically perfect flatness, but relies on **piecewise "stitching"** of functions. This creates impulsive "kicks" (discontinuities) in higher-order derivatives (Jerk/Snap) and makes hardware optimization difficult.

### 🔬 Theoretical Basis: CMST Theory
This window is an implementation of **CMST (Cosh Moment Sturm Transform)**. Unlike standard windows which are often heuristic curve-fits, the CMST window is constructed as a **Geometric Mollifier** with three rigorous guarantees:

1.  **Analytically Sound ($C^\infty$):**
    The function belongs to the Gevrey class of regularity ($s=2$), ensuring super-algebraic decay in the frequency domain. It is infinitely differentiable ($C^\infty$) with no discontinuities in any derivative. This eliminates the "spectral ringing" and mechanical jerk caused by piecewise functions like the Planck-taper or Tukey window.

2.  **Zero-Preserving :**
    Derived from CMST theory, the kernel guarantees the preservation of realness in the signal chain. It does not introduce artificial complex roots (phantom oscillations) into the passband.

3.  **Interlace-Preserving Transform:**
    For all derivatives, the window acts as a variation-diminishing operator. It preserves the root-interlacing structure of the underlying signal, ensuring that derivative noise is bounded and geometric topology is maintained even at the boundaries.

### 💡 The Formula (CMST)
We utilize a compensated log-concave mollifier that cancels low-order curvature to achieve Flatness:

$$
w(t) = \exp\left(1+t^p - \frac{1}{1-t^p}\right), 
$$

where p is even.

* **Compensating Term $t^p$:** Cancels the Gaussian curvature at the origin, extending the "Table-Top" flatness to order $2n$.
* **Mollifier Term $(-1/(1-t^p))$ :** Enforces strict compact support with essential singularities at the boundaries, ensuring all derivatives decay to zero smoothly.


### 🎛️ Tunable Flatness (p-Control)

Unlike traditional windows which are locked to a single profile (e.g., Hann, Blackman), the CMST window is a parametric family. The power parameter (p) allows you to tune the window's behavior to match your specific engineering constraint:

    Mode A: The "Brick Wall" (p=6 or higher)
        Goal: Maximal Amplitude Accuracy.

        Behavior: The window remains effectively flat (>0.99) for over 70% of the duration, ensuring that signals are not attenuated in the center.

    Mode B: The "Silencer" (p=2)
        Goal: Maximal Spectral Purity.

        Behavior: The window converges to an analytically smooth Gaussian-like profile. This sacrifices the "flat top" to achieve significantly faster side-lobe decay, often over 100dB improvement, diving into the noise floor deeper than standard piecewise functions like the Planck-taper.
        Generally I have found p=2 gives the best results.

### 📉 Tunable Flatness
![Tunable shape](tunable_flatness_time_only.png)

### 📉 Spectral Leakage Comparison
![Comparison](VsOtherFilters.png)

Note on Precision Limits: The CMST response (Blue) was calculated using 100 digit precision arithmetic to demonstrate the asymptotic behavior beyond standard 64-bit machine limits (~ -320 dB).
Polynomial Windows (e.g., Blackman-Harris, Planck): Their leakage floors are theoretical. Even with infinite precision, they would not drop much further and the graph wouldnt change much


### 📉 Performance Analysis: 

**Detection of a -100 dB weak signal (1.3 kHz) next to a strong carrier (1.0 kHz).**

The standard Planck-taper (Red) buries the target in spectral leakage. The CMST window (Blue) resolves it clearly with >30 dB of headroom
![Discrimination](weak_signal_comparison.png)

**Transform of a Sinc function**
![Sinc Transform](spectral_comparison_db.png)

Lets remember that when I say these should be a box, we are using a log scale, visually they are all boxes.

Total Integrated Leakage (0.05 - 0.07 Hz):

Planck Taper:  -20.86 dB

7-term BH:     -16.21 dB

CMST (p=2):    -28.44 dB

Improvement:   7.58 dB

**Planet sim.**
![Planet sim](planets.png)

With the planet dB at -80 dB and a noise floor of -90 dB.  Slightly contrived, but the windows I am comparing against are good!
Note that not only has Planck lost 4 planets, but it has also created non existant mountains (side lobes).  We are looking at this in crazy detail, if we were back in the real world, the -80dB planets would not show up on your screen, your eyes couldnt see them.  Those "mountains" are invisible as well. 80dB is the ratio between the Empire State Building and a golf ball.


### 📉 Resolution Law: 
As part of this work we produce a resolution law for CMST where for resolution of two signals in terms of bins, namely

$$m > \left\lceil \frac{(\ln R)^2}{\pi} \right\rceil$$

Where:
* **$R$** is the linear ratio of amplitudes (e.g., for -100 dB, $R = 10^5$).
* **$m$** is the distance between the signals in bins.

![Resolution Law](resolution_law.png)


### 📉 The Math: 
Behind all of this there is a CMST theory paper here [CMST](cmst.pdf)


### 📦 Installation

```bash
git clone https://github.com/aronp/CMST.git

cd CMST

pip install .
