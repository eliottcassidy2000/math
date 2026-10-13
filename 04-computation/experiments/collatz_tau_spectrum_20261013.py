"""
collatz_tau_spectrum_20261013.py -- the log-periodic bias of the parity of the total stopping time.
Bin n by 1/16 octave from 2^18 to 2^29 (176 bins), take the mean of (-1)^tau per bin, and estimate the spectrum in
x = log_2 n (cycles per octave) by least squares on a frequency grid; report the dominant frequencies, amplitudes, and
the per-octave rms amplitude (does the bias decay?).  Also the logarithmic mean (sum t(n)/n / log N) as a control.
Run: python3 collatz_tau_spectrum_20261013.py
"""
import time
import numpy as np
from numba import njit

@njit(cache=True)
def tau_array16(N):
    tau = np.zeros(N + 8, dtype=np.int16)
    for n in range(2, N + 8):
        m = n; steps = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]
    return tau

t0 = time.time()
N = 1 << 29
tau = tau_array16(N)
t = np.where(tau % 2 == 1, -1.0, 1.0); t[0] = 0.0
print(f'tau to 2^29 in {time.time()-t0:.1f}s')
J0, J1, R = 18, 29, 16
edges = [int(round(2 ** (J0 + k / R))) for k in range((J1 - J0) * R + 1)]
xs = np.array([(np.log2(edges[k]) + np.log2(edges[k + 1])) / 2 for k in range(len(edges) - 1)])
ms = np.array([t[edges[k]:edges[k + 1]].mean() for k in range(len(edges) - 1)])
print(f'{len(ms)} bins of 1/{R} octave; overall mean {ms.mean():+.5f}; rms {np.sqrt((ms**2).mean()):.5f}')
print('per-octave rms of the binned bias:')
for j in range(J0, J1):
    sel = (xs >= j) & (xs < j + 1)
    print(f'  [2^{j}, 2^{j+1}): rms {np.sqrt((ms[sel]**2).mean()):.5f}  mean {ms[sel].mean():+.5f}')
# least-squares spectrum on a frequency grid (cycles per octave)
freqs = np.linspace(0.05, 4.0, 791)
power = []
for f in freqs:
    A = np.column_stack([np.cos(2 * np.pi * f * xs), np.sin(2 * np.pi * f * xs)])
    coef, *_ = np.linalg.lstsq(A, ms, rcond=None)
    power.append(np.hypot(*coef))
power = np.array(power)
order = np.argsort(-power)
peaks = []
for i in order:
    f = freqs[i]
    if all(abs(f - g) > 0.08 for g, _ in peaks):
        peaks.append((f, power[i]))
    if len(peaks) >= 8:
        break
print('dominant frequencies (cycles per octave) and amplitudes:')
for f, a in peaks:
    print(f'  f = {f:.3f}  period {1/f:.3f} octaves  amplitude {a:.5f}   f*log2(3) mod 1 = {(f*np.log2(3)) % 1:.3f}')
# joint fit of the top 4 and the residual
A = np.column_stack([fn(2 * np.pi * f * xs) for f, _ in peaks[:4] for fn in (np.cos, np.sin)] + [np.ones_like(xs)])
coef, *_ = np.linalg.lstsq(A, ms, rcond=None)
res = ms - A @ coef
print(f'joint fit of top 4 modes: residual rms {np.sqrt((res**2).mean()):.5f} vs data rms {np.sqrt((ms**2).mean()):.5f}; constant term {coef[-1]:+.5f}')
# logarithmic mean control
w = 1.0 / np.arange(1, N + 1)
for J in (20, 24, 29):
    M = 1 << J
    print(f'logarithmic mean to 2^{J}: {float((t[1:M+1] * w[:M]).sum() / np.log(M)):+.5f}')
print('done')
