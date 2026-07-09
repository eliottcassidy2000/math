"""
lrc14_L2_route_test_opus_S172.py   (opus-2026-07-09-S172)

Does the L2 route rescue the a-priori residual bound the TV/L1 route could not?  Cauchy-Schwarz with
weight m^2 on the resonant sum R = sum_{Vmax|m} What_1d(m):
   |R| <= sum_{n>=1} 2|What_1d(nVmax)|
        <= 2 (sum_{n>=1} |nVmax*What_1d(nVmax)|^2)^{1/2} (sum_{n>=1} 1/(nVmax)^2)^{1/2}
        <= 2 (E[(W')^2])^{1/2} * (zeta(2))^{1/2} / Vmax          [Parseval: sum|m What(m)|^2 = ||W'||_2^2]
So the L2 a-priori bound is  B_L2 = 2*sqrt(E[(W')^2]*zeta(2))/Vmax.  It CLOSES (B_L2<main) iff
E[(W')^2] is small enough -- the over-covering (k/7>1) makes W'=0 most of the time, so MAYBE E[(W')^2]
<< spread^2 and the L2 route works where TV (~spread^2) failed.  DECISIVE TEST.
"""
import numpy as np
from math import gcd, pi

K = 13
MAIN = (6.0 / 7.0) ** K
Z2 = pi ** 2 / 6


def W_series(E, x):
    Earr = np.asarray(E, float)
    Ph = np.mod(np.outer(x, Earr), 1.0); Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:, -1] + Ph[:, 0])[:, None]], axis=1)
    return np.clip(g - 1.0 / 7, 0, None).sum(axis=1)


def E_Wprime2(E, NG):
    """E[(W')^2] via forward difference on a fine grid."""
    x = np.arange(NG) / NG
    W = W_series(E, x)
    Wp = (np.roll(W, -1) - W) * NG        # W' ~ dW/dx
    return float((Wp ** 2).mean())


def R_abs_conv(E, Vmax, Lmax=250):
    N = Lmax * Vmax
    x = np.arange(N) / N
    F = np.fft.rfft(W_series(E, x)) / N
    ra = 0.0; rs = 0.0
    for n in range(1, Lmax):
        m = n * Vmax
        if m >= F.shape[0]:
            break
        ra += 2 * abs(F[m]); rs += 2 * F[m].real
    return ra, rs


def prim(E):
    d = 0
    for i in range(len(E) - 1):
        d = gcd(d, E[i + 1] - E[i])
    return d == 1


import random
rng = random.Random(1724)
print("=" * 96)
print(f"L2 route test: B_L2 = 2*sqrt(E[(W')^2]*zeta2)/Vmax  vs main=(6/7)^{K}={MAIN:.5f}")
print("=" * 96)
print(f"  {'spread':>7} {'Vmax':>6} {'E_Wp2':>10} {'sqrt/spr':>9} {'B_L2':>9} {'B_L2/main':>10} "
      f"{'|R|actual':>10} {'B_L2<main?':>10}")
for s in [30, 50, 80, 120, 180, 260, 380]:
    bestB = 0; row = None
    for _ in range(6):
        mid = sorted(rng.sample(range(1, s), K - 2)); E = [0] + mid + [s]
        if len(set(E)) != K or not prim(E):
            continue
        NG = min(60000, 500 * s)
        ew2 = E_Wprime2(E, NG)
        Vmax = s + 1
        B = 2 * np.sqrt(ew2 * Z2) / Vmax
        if B > bestB:
            bestB = B; row = (Vmax, ew2, B, E)
    Vmax, ew2, B, E = row
    ra, rs = R_abs_conv(E, Vmax)
    print(f"  {s:>7} {Vmax:>6} {ew2:>10.2f} {np.sqrt(ew2)/s:>9.4f} {B:>9.4f} {B/MAIN:>10.2f} "
          f"{abs(rs):>10.5f} {'YES' if B < MAIN else 'NO':>10}")
print()
print("  READING: if E[(W')^2]~spread^2 (sqrt/spr ~ const) then B_L2 ~ const*sqrt(zeta2) -> NOT <main")
print("  (L2 route ALSO hits the wall => |R|<main is ARITHMETIC cancellation, not analytic magnitude).")
print("  if E[(W')^2]<<spread^2 (over-covering kills W') then B_L2->0 and the L2 route CLOSES it.")
print("=" * 96)
