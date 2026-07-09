"""
lrc14_Rabs_crossover_opus_S172.py   (opus-2026-07-09-S172)

HONEST reconciliation with kps-S98 (and correction of my own S171 over-claim).  kps-S98 finds the
absolute resonant sum R_abs is NOT uniformly < main=(6/7)^k: ~0.6 large spread, but >1 at small-moderate
spread (1.55@s50, 1.12@s70), while signed |R|~0.1 stays small => small-spread smallness is CANCELLATION
(the Mertens wall), not an absolute bound.  My S171 adversarial reported max 0.40 -- likely an L=64 FFT
TRUNCATION artifact.  This recomputes R_abs (per-frequency, CONVERGED) vs spread to locate the crossover
R_abs=main, confirm large-spread R_abs<main (a-priori-closable) vs small-spread R_abs>main (Mertens/
exhaustion).  R_abs = sum_{n>=1} 2|What_1d(nVmax)|, What_1d(m)=1-D Fourier coeff of W(x).
"""
import numpy as np
from math import gcd
from functools import reduce

K = 13
MAIN = (6.0 / 7.0) ** K


def W_series(E, x):
    Earr = np.asarray(E, float)
    Ph = np.mod(np.outer(x, Earr), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return np.clip(g - 1.0 / 7, 0, None).sum(axis=1)


def R_abs_converged(E, Vmax, Lmax=400):
    """R_abs = sum_{n>=1} 2|What_1d(nVmax)|, using a fine grid N=Lmax*Vmax; also |R| signed and #terms."""
    N = Lmax * Vmax
    x = np.arange(N) / N
    F = np.fft.rfft(W_series(E, x)) / N
    Rabs = 0.0; Rsig = 0.0; nz = 0
    for n in range(1, Lmax):
        m = n * Vmax
        if m >= F.shape[0]:
            break
        a = 2 * abs(F[m])
        Rabs += a; Rsig += 2 * F[m].real
        if a > 1e-6:
            nz = n
    return Rabs, Rsig, nz


def prim(E):
    d = 0
    for i in range(len(E) - 1):
        d = gcd(d, E[i + 1] - E[i])
    return d == 1


import random
rng = random.Random(172)
print("=" * 92)
print(f"R_abs (per-frequency, CONVERGED) vs spread.  main=(6/7)^{K}={MAIN:.5f}  (dissociated k=13)")
print("=" * 92)
print(f"  {'spread':>7} {'Vmax':>6} {'R_abs':>9} {'R_abs/main':>10} {'|R_sig|':>9} {'|Rs|/main':>9} {'ntail':>6} {'<main?':>7}")
for s in [30, 40, 50, 60, 70, 90, 120, 160, 220, 300, 400]:
    # sample a few dissociated clusters at this spread, report the WORST R_abs
    worst = 0.0; row = None
    for _ in range(8):
        mid = sorted(rng.sample(range(1, s), K - 2))
        E = [0] + mid + [s]
        if len(set(E)) != K or not prim(E):
            continue
        Vmax = s + 1
        ra, rs, nt = R_abs_converged(E, Vmax, Lmax=300)
        if ra > worst:
            worst = ra; row = (Vmax, ra, rs, nt)
    if row is None:
        continue
    Vmax, ra, rs, nt = row
    print(f"  {s:>7} {Vmax:>6} {ra:>9.5f} {ra/MAIN:>10.3f} {abs(rs):>9.5f} {abs(rs)/MAIN:>9.3f} {nt:>6} "
          f"{'YES' if ra < MAIN else 'NO':>7}")
print()
print("  READING: if R_abs/main >1 at small-moderate spread and <1 (=>~0.6) at large spread, my S171")
print("  'R_abs<main uniformly, no cancellation' was an OVER-CLAIM (truncation) -- converges with kps-S98:")
print("  small-spread dissociated is the Mertens wall (signed |R| small by cancellation, R_abs>main),")
print("  closed-form absolute bound only LARGE-spread; small-spread => LEM-013 exhaustion.")
print("=" * 92)
