"""
lrc14_Egrid_residual_absolute_opus_S171.py   (opus-2026-07-09-S171)

DECISIVE diagnostic for kps-S96's E_grid[W]>0 existence route (LRCEgridExistence.lean).  That route
closes the good period a-priori IF |R| < (6/7)^k, R = sum_{Vmax|m, m!=0} What(m) the resonance residual
(What = Fourier coeff of the uncovered measure W(x)=sum_i (gap_i-1/7)_+).  kps calls it "Mertens-safe (a
count, not a cancellation)".  THE QUESTION this settles: is R controlled by an ABSOLUTE bound
  R_abs := sum_{n>=1} 2|What(nVmax)|  <  (6/7)^k        [=> a-priori via decay+count, NO cancellation]
or only by SIGNED cancellation R_signed := sum_{n>=1} 2 Re What(nVmax), with R_abs >> |R_signed|
[=> Mertens-hard]?  W is CONTINUOUS piecewise-linear (unlike the sharp indicator), so What(m)~1/m^2
(opus-S170 alpha=2), which should make R_abs itself small -- the smooth route makes kps's route ABSOLUTE.

Computes, for dissociated and 7-STRUCTURED (diffs=0 mod 7, mac-mini MISTAKE-128) k=13 clusters at the
critical Vmax: (6/7)^13, R_signed (== E_grid[W]-(6/7)^13, exact ruler), R_abs (FFT), the ratio
R_abs/|R_signed| (=1 => absolute, >>1 => cancellation), R_abs/(6/7)^13 (<1 => absolute bound closes it),
and the What decay exponent alpha.
"""
import numpy as np
from math import gcd
from functools import reduce

INV7 = 1.0 / 7.0
K = 13
MAIN = (6.0 / 7.0) ** K


def W_series(E, x):
    """uncovered measure W(x)=sum_i (gap_i-1/7)_+ for array x."""
    E = np.asarray(E, float)
    Ph = np.mod(np.outer(x, E), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return np.clip(g - INV7, 0, None).sum(axis=1)


def Egrid(E, Vmax):
    j = np.arange(Vmax); x = j / Vmax
    return float(W_series(E, x).mean())


def residuals(E, Vmax, L=96):
    """R_signed (exact ruler) and R_abs (FFT over m=nVmax) and What decay alpha."""
    R_signed = Egrid(E, Vmax) - MAIN
    N = L * Vmax
    x = np.arange(N) / N                      # grid j/N (matches ruler j/Vmax at n*Vmax)
    f = W_series(E, x)
    F = np.fft.rfft(f) / N                     # What(m) ~ F[m]
    R_abs = 0.0; R_sig_fft = 0.0
    for n in range(1, L):
        m = n * Vmax
        if m < F.shape[0]:
            R_abs += 2 * abs(F[m])
            R_sig_fft += 2 * F[m].real
    # decay alpha of |What(m)| over a mid band
    ms = np.arange(2, min(4000, F.shape[0]))
    lF = np.log(np.maximum(np.abs(F[ms]), 1e-15)); lm = np.log(ms)
    alpha = -np.polyfit(lm, lF, 1)[0]
    return R_signed, R_abs, R_sig_fft, alpha


def longest_ap(E):
    S = set(E); E = sorted(E); best = 2
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            d = E[j] - E[i]; L = 2; nx = E[j] + d
            while nx in S: L += 1; nx += d
            best = max(best, L)
    return best


print("=" * 100)
print("E_grid residual: is |R| ABSOLUTE (R_abs<(6/7)^k, a-priori) or CANCELLATION (R_abs>>|R_signed|)?")
print(f"  k={K}, main term (6/7)^{K} = {MAIN:.5f}")
print("=" * 100)
families = {
    "dissociated Sidon      ": [1, 2, 5, 11, 22, 33, 45, 60, 78, 95, 110, 130, 150],
    "dissociated random     ": sorted(np.random.RandomState(3).choice(range(1, 120), 13, replace=False).tolist()),
    "7-structured (=1 mod 7) ": [1 + 7 * i for i in range(13)],
    "7-structured (=3 mod 7) ": [3 + 7 * i for i in range(13)],
    "near-AP 7*{1..10}+1,2,3 ": sorted(set([7 * i for i in range(1, 11)] + [1, 2, 3])),
}
print(f"  {'family':>26} {'L':>3} {'Vmax':>6} {'R_signed':>9} {'R_abs':>8} {'R_abs/|Rs|':>10} "
      f"{'R_abs/main':>10} {'alpha':>6} {'|Rs|<main?':>10} {'R_abs<main?':>11}")
for name, E in families.items():
    E = sorted(set(E))
    if len(E) != 13:
        continue
    spread = max(E)
    Vmax = spread + 1                       # binding end of the critical window [s+1, 7s/6]
    Rs, Ra, Rsf, alpha = residuals(E, Vmax)
    L = longest_ap(E)
    ratio = Ra / max(abs(Rs), 1e-9)
    print(f"  {name:>26} {L:>3} {Vmax:>6} {Rs:>9.5f} {Ra:>8.5f} {ratio:>10.1f} "
          f"{Ra / MAIN:>10.3f} {alpha:>6.2f} {('YES' if abs(Rs) < MAIN else 'NO'):>10} "
          f"{('YES' if Ra < MAIN else 'NO'):>11}")

print()
print("  READING: R_abs/|Rs| ~ 1 => the residual is ABSOLUTELY bounded (no cancellation) => kps's")
print("  |R|<(6/7)^k is a-priori-provable via [What~1/m^alpha decay]x[resonance count]. R_abs>>|Rs|")
print("  => Mertens cancellation (hard). R_abs<main (last col) => the ABSOLUTE bound alone closes it.")
print("=" * 100)
