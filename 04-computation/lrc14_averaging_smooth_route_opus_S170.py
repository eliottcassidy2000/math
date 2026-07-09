"""
lrc14_averaging_smooth_route_opus_S170.py   (opus-2026-07-09-S170)

THE SMOOTH AVERAGING ROUTE to good-period existence -- sidesteps the Mertens/L2 signed cancellation
of the sharp-indicator Erdos-Turan sum (klein MISTAKE-127 / opus-S170 identity script).

KEY: a lonely (good) period is j with maxgap_13{frac(e_i j/Vmax)} > 1/7 (observer in the max gap of
all 13 runners).  The AVERAGING route (kps-S95): if the ruler-grid MEAN E_j[maxgap] > 1/7 then some j
is good (max >= mean).  Crucially maxgap(x) is CONTINUOUS piecewise-linear (circular gaps vary
continuously), so its Fourier coefficients decay like 1/m^2 (kink, not jump) -- UNLIKE the sharp
good-set indicator (1/m, L1-divergent, needs cancellation).  Hence the averaging discrepancy
  E_j[maxgap] - E_x[maxgap] = sum_{n!=0} maxgap^(nVmax) (-1)^n,   |.| <= C/Vmax^2  (CONVERGES),
so E_j[maxgap] -> E_x[maxgap] a-priori, and existence reduces to the DILATION-INVARIANT inequality
  E_x[maxgap] > 1/7  (+ margin),
with NO cancellation.  DECISIVE TEST: does E_x[maxgap] > 1/7 hold for the EXTREMAL families
(tight AP {1..13}, dilated AP, near-AP dilated-AP+outliers), the ones that break arc-count?

This script: (A) E_x[maxgap_13] (continuous, dilation-invariant) and margin over 1/7 for extremal +
dissociated 13-runner families; (B) the maxgap Fourier decay exponent (~2 => convergent resonant
sum); (C) the averaging discrepancy E_j - E_x vs Vmax (~1/Vmax^2), => a-priori E_j[maxgap]>1/7 for
Vmax >= V0.
"""
import numpy as np

INV7 = 1.0 / 7.0


def maxgap_series(E, x):
    """maxgap_{|E|}{frac(e_i x)} for array x (circular largest gap)."""
    E = np.asarray(E, float)
    Ph = np.mod(np.outer(x, E), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return g.max(axis=1)


def Ex_maxgap(E, N=2_000_000):
    x = (np.arange(N) + 0.5) / N
    return float(maxgap_series(E, x).mean())


def Ej_maxgap(E, Vmax):
    j = np.arange(Vmax); x = (j + 0.5) / Vmax
    return float(maxgap_series(E, x).mean())


def fourier_decay(E, N=1 << 20):
    """fit |maxgap^(m)| ~ m^-alpha over a mid band; return alpha."""
    x = (np.arange(N) + 0.5) / N
    f = maxgap_series(E, x)
    F = np.abs(np.fft.rfft(f - f.mean()) / N)
    ms = np.arange(1, F.shape[0])
    # mid band 8..2000, average |F| in dyadic bins to fit slope
    lo, hi = 8, 2000
    band = (ms >= lo) & (ms <= hi)
    lm = np.log(ms[band]); lF = np.log(np.maximum(F[1:][band], 1e-15))
    A = np.polyfit(lm, lF, 1)
    return -A[0]


print("=" * 92)
print("(A) E_x[maxgap_13] (continuous, dilation-invariant) and margin over 1/7 = 0.142857")
print("=" * 92)
families = {
    "tight AP {1..13}          ": list(range(1, 14)),
    "dilated AP 5*{1..13}      ": [5 * i for i in range(1, 14)],
    "near-AP 7*{1..10}+{1,2,3} ": sorted(set([7 * i for i in range(1, 11)] + [1, 2, 3])),
    "near-AP 11*{1..10}+{1,2,3}": sorted(set([11 * i for i in range(1, 11)] + [1, 2, 3])),
    "dissociated Sidon-ish     ": [1, 2, 5, 11, 22, 33, 45, 60, 78, 95, 110, 130, 150],
    "random primitive          ": sorted(np.random.RandomState(7).choice(range(1, 90), 13, replace=False)),
}
print(f"  {'family':>28} {'k':>3} {'E_x[maxgap]':>12} {'ratio/(1/7)':>12} {'>1/7?':>6}")
Ex_store = {}
for name, E in families.items():
    ex = Ex_maxgap(E)
    Ex_store[name] = (E, ex)
    print(f"  {name:>28} {len(E):>3} {ex:>12.5f} {ex/INV7:>12.4f} {('YES' if ex > INV7 else 'NO'):>6}")

print()
print("=" * 92)
print("(B) maxgap Fourier decay exponent alpha (|maxgap^(m)| ~ m^-alpha); alpha~2 => 1/m^2, convergent")
print("=" * 92)
for name in ["tight AP {1..13}          ", "near-AP 7*{1..10}+{1,2,3} ", "dissociated Sidon-ish     "]:
    E = families[name]
    a = fourier_decay(E)
    print(f"  {name:>28}: alpha = {a:.3f}")

print()
print("=" * 92)
print("(C) averaging discrepancy E_j[maxgap]-E_x[maxgap] vs Vmax (near-AP d*{1..10}+{1,2,3})")
print("=" * 92)
print(f"  {'d':>4} {'Vmax':>6} {'E_x':>9} {'E_j':>9} {'|E_j-E_x|':>10} {'*Vmax^2':>9} {'E_j>1/7?':>8}")
for d in (5, 10, 20, 40, 80, 160):
    E = sorted(set([d * i for i in range(1, 11)] + [1, 2, 3]))
    Vmax = max(E) + 1
    ex = Ex_maxgap(E); ej = Ej_maxgap(E, Vmax)
    disc = abs(ej - ex)
    print(f"  {d:>4} {Vmax:>6} {ex:>9.5f} {ej:>9.5f} {disc:>10.6f} {disc*Vmax**2:>9.2f} "
          f"{('YES' if ej > INV7 else 'NO'):>8}")

print()
print("  READING: if (A) E_x[maxgap]>1/7 with a UNIFORM margin (incl. the extremal AP families that")
print("  break arc-count) and (B) alpha~2 (maxgap Fourier ~1/m^2), then (C) |E_j-E_x|=O(1/Vmax^2)->0,")
print("  so E_j[maxgap]>1/7 a-priori for Vmax>=V0 => a GOOD PERIOD EXISTS -- NO cancellation, NO arc-count.")
print("=" * 92)
