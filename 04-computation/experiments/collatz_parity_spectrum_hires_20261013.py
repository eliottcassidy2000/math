"""
collatz_parity_spectrum_hires_20261013.py -- per-octave high-resolution spectrum of the parity of the total stopping time.
tau to N = 2^NL (int16); bins of 1/R octave; for each octave j, least-squares amplitude of the binned parity mean at a
frequency grid (cycles per octave) up to the Nyquist R/2; peaks tracked across octaves (amplitude, phase); then a joint
damped-sinusoid fit over octaves J0..J1 of the dominant low mode with FREE frequency, compared with the two candidates
f = 1/(2 - log2 3) = 2.4095 (half drift) and f = 5/2 (signed resonance at 5 log2 3 ~ 8).
Run: python3 collatz_parity_spectrum_hires_20261013.py <NL> <R>
"""
import sys, time
import numpy as np
from numba import njit

NL = int(sys.argv[1]) if len(sys.argv) > 1 else 30
R = int(sys.argv[2]) if len(sys.argv) > 2 else 64
THETA = np.log2(3.0)

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

@njit(cache=True)
def bin_means(tau, J0, J1, R):
    nb = (J1 - J0) * R
    sums = np.zeros(nb); cnts = np.zeros(nb)
    lo = 1 << J0; hi = 1 << J1
    for n in range(lo, hi):
        x = np.log2(n) - J0
        k = int(x * R)
        if k >= nb:
            k = nb - 1
        v = -1.0 if (tau[n] & 1) else 1.0
        sums[k] += v; cnts[k] += 1
    return sums / np.maximum(cnts, 1), cnts

t0 = time.time()
N = 1 << NL
tau = tau_array16(N)
print(f'tau to 2^{NL} in {time.time()-t0:.1f}s; max {int(tau.max())}', flush=True)
J0, J1 = 16, NL
ms, cnts = bin_means(tau, J0, J1, R)
xs = J0 + (np.arange(len(ms)) + 0.5) / R
print(f'{len(ms)} bins of 1/{R} octave from 2^{J0} to 2^{J1}; min count per bin {int(cnts.min())}')

def ls_amp(x, y, f):
    A = np.column_stack([np.cos(2 * np.pi * f * x), np.sin(2 * np.pi * f * x), np.ones_like(x)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return np.hypot(c[0], c[1]), np.arctan2(-c[1], c[0]), c[2]

freqs = np.arange(0.25, R / 2, 0.05)
print('\n# per-octave spectra: top 6 peaks (f, amplitude), and the amplitudes at the candidate frequencies')
cands = [1 / (2 - THETA), 2.5, 6 / (2 - THETA), 14.5, 11 / (2 - THETA), 26.5]
print('   candidates:', ' '.join(f'{c:.4f}' for c in cands))
rows = []
for j in range(J0, J1):
    sel = (xs >= j) & (xs < j + 1)
    x = xs[sel]; y = ms[sel]
    amps = np.array([ls_amp(x, y, f)[0] for f in freqs])
    order = np.argsort(-amps)
    peaks = []
    for i in order:
        if all(abs(freqs[i] - g) > 0.3 for g, _ in peaks):
            peaks.append((freqs[i], amps[i]))
        if len(peaks) >= 6:
            break
    ca = [ls_amp(x, y, c)[0] for c in cands]
    rows.append((j, peaks, ca))
    print(f' 2^{j:2d}: ' + ' '.join(f'({f:.2f},{a:.4f})' for f, a in peaks) + ' | cand ' + ' '.join(f'{a:.4f}' for a in ca), flush=True)

print('\n# joint fit over octaves 20..J1 of the low mode: damped sinusoid A exp(-g x) cos(2 pi f x + phi), free f near 2.45')
sel = xs >= 20
x = xs[sel]; y = ms[sel]
def resid(f, g):
    w = np.exp(-g * x)
    A = np.column_stack([w * np.cos(2 * np.pi * f * x), w * np.sin(2 * np.pi * f * x), np.ones_like(x)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c
    return float((r ** 2).sum()), c
best = None
for f in np.arange(2.30, 2.60, 0.0025):
    for g in np.arange(-0.02, 0.20, 0.005):
        s, c = resid(f, g)
        if best is None or s < best[0]:
            best = (s, f, g, c)
s0 = float((y ** 2).sum())
print(f'best: f = {best[1]:.4f}, decay g = {best[2]:.3f} per octave (factor {np.exp(-best[2]):.3f}), residual/total = {best[0]/s0:.4f}')
for f in (1 / (2 - THETA), 2.5):
    bs = None
    for g in np.arange(-0.02, 0.20, 0.0025):
        s, c = resid(f, g)
        if bs is None or s < bs[0]:
            bs = (s, g)
    print(f'  fixed f = {f:.4f}: best decay g = {bs[1]:.4f} (factor {np.exp(-bs[1]):.3f}), residual/total = {bs[0]/s0:.4f}')
# the same for the second mode near 14.5
best2 = None
for f in np.arange(14.30, 14.70, 0.0025):
    for g in np.arange(-0.05, 0.20, 0.005):
        s, c = resid(f, g)
        if best2 is None or s < best2[0]:
            best2 = (s, f, g)
print(f'second mode: best f = {best2[1]:.4f}, decay g = {best2[2]:.3f} (factor {np.exp(-best2[2]):.3f}); candidates {6/(2-THETA):.4f} and 14.5')
best3 = None
for f in np.arange(26.30, 26.70, 0.0025):
    for g in np.arange(-0.05, 0.20, 0.005):
        s, c = resid(f, g)
        if best3 is None or s < best3[0]:
            best3 = (s, f, g)
print(f'third mode: best f = {best3[1]:.4f}, decay g = {best3[2]:.3f} (factor {np.exp(-best3[2]):.3f}); candidates {11/(2-THETA):.4f} and 26.5')
np.save('collatz_parity_spectrum_hires_20261013_bins.npy', np.vstack([xs, ms, cnts]))
print('done')
