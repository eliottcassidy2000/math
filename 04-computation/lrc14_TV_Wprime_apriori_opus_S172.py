"""
lrc14_TV_Wprime_apriori_opus_S172.py   (opus-2026-07-09-S172)

The decisive number for the OCF/decay TRANSPORT.  The good-period E_grid residual satisfies
  |R| <= R_abs^{per-freq} = sum_{n>=1} 2|What_1d(nVmax)|,   What_1d = 1-D Fourier coeff of W(x).
Since W = sum(gap-1/7)_+ is CONTINUOUS piecewise-linear (W' is BV), the rigorous a-priori bound is
  |What_1d(m)| <= TV(W')/(2pi m)^2   =>   R_abs <= TV(W')/(12 Vmax^2).
So the transport gives an a-priori R_abs < main=(6/7)^k  IFF  TV(W') < 12*main*Vmax^2 = 1.62*Vmax^2.

CRUX: how does TV(W') scale?  Naive worst case ~56*spread^2 (all C(13,2) pairs) would give R_abs<=4.7
(useless).  BUT only the <=6 gaps EXCEEDING 1/7 contribute to W' -- so TV(W') may be O(spread), giving
R_abs <= O(1/spread) -> 0 (closes at large spread).  This computes TV(W') (via |second difference|),
its scaling exponent, the a-priori bound TV/(12Vmax^2), and compares to the ACTUAL R_abs and main.
"""
import numpy as np
from math import gcd

K = 13
MAIN = (6.0 / 7.0) ** K


def W_series(E, x):
    Earr = np.asarray(E, float)
    Ph = np.mod(np.outer(x, Earr), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return np.clip(g - 1.0 / 7, 0, None).sum(axis=1)


def TV_Wprime(E, NG):
    """TV(W') = integral |W''| ~ sum |W(x+h)-2W(x)+W(x-h)|/h over a fine periodic grid."""
    x = np.arange(NG) / NG
    W = W_series(E, x)
    Wpp = np.abs(np.roll(W, -1) - 2 * W + np.roll(W, 1))       # * (1/h^2), times h = *NG
    return float(Wpp.sum() * NG)     # sum |second diff| / h,  h=1/NG => * NG


def R_abs_conv(E, Vmax, Lmax=250):
    N = Lmax * Vmax
    x = np.arange(N) / N
    F = np.fft.rfft(W_series(E, x)) / N
    ra = 0.0
    for n in range(1, Lmax):
        m = n * Vmax
        if m >= F.shape[0]:
            break
        ra += 2 * abs(F[m])
    return ra


def prim(E):
    d = 0
    for i in range(len(E) - 1):
        d = gcd(d, E[i + 1] - E[i])
    return d == 1


import random
rng = random.Random(1723)
print("=" * 100)
print(f"TV(W') scaling + the a-priori transport bound.  main=(6/7)^{K}={MAIN:.5f}, need TV<1.62*Vmax^2")
print("=" * 100)
print(f"  {'spread':>7} {'Vmax':>6} {'TVWp':>9} {'TV/spr':>8} {'TV/spr^2':>9} {'TV/(12V^2)':>11} "
      f"{'R_abs':>8} {'bound<main?':>11} {'Rabs<main?':>10}")
sprs = [30, 50, 80, 120, 180, 260, 380]
tv_list = []; sp_list = []
for s in sprs:
    # worst (max TV and max R_abs) over a few dissociated clusters
    bestTV = 0; bestRA = 0; E_at = None
    for _ in range(6):
        mid = sorted(rng.sample(range(1, s), K - 2)); E = [0] + mid + [s]
        if len(set(E)) != K or not prim(E):
            continue
        NG = min(40000, 400 * s)
        tv = TV_Wprime(E, NG)
        if tv > bestTV:
            bestTV = tv; E_at = E
    ra = R_abs_conv(E_at, s + 1)
    Vmax = s + 1
    bound = bestTV / (12 * Vmax ** 2)
    tv_list.append(bestTV); sp_list.append(s)
    print(f"  {s:>7} {Vmax:>6} {bestTV:>9.2f} {bestTV/s:>8.3f} {bestTV/s**2:>9.4f} {bound:>11.4f} "
          f"{ra:>8.4f} {'YES' if bound < MAIN else 'NO':>11} {'YES' if ra < MAIN else 'NO':>10}")
# scaling exponent
if len(tv_list) >= 2:
    a = np.polyfit(np.log(sp_list), np.log(tv_list), 1)[0]
    print(f"\n  FIT: TV(W') ~ spread^{a:.2f}  (exponent <2 => TV/(12Vmax^2) -> 0 as spread grows => a-priori)")
print()
print("  READING: if TV(W')~spread^1 (only the few >1/7 gaps contribute) then TV/(12Vmax^2)~1/(12 spread)")
print("  -> 0, and the TV bound gives R_abs<main a-priori for large spread. If TV~spread^2, the crude TV")
print("  bound is useless (~const) and R_abs<main is genuine CANCELLATION (kps-S98 Mertens wall).")
print("=" * 100)
