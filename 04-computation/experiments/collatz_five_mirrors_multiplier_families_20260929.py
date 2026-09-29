#!/usr/bin/env python3
"""Does the maximal primitive Fourier coefficient stay on the pure powers of two?  (S20 addendum.)

For every unit multiplier u the family t = u 2^j (j in Z) is closed under the exact frequency recursion
(u 2^j -> u 2^(j-a)), so m_n^u(j) := mu_hat_n(u 2^j mod 3^n) can be followed level by level as in
collatz_five_mirrors_powers_of_two_20260929.py.  We compare max_j |m_h^u(j)| for the odd multipliers u <= UMAX
coprime to 3 (u = 1 is the pure family; |mu_hat(-t)| = |mu_hat(t)| so u > 0 suffices) at levels h <= HMAX.
For h <= 18 every unit is u 2^j for some u <= ..., so the FFT maxima are reproduced only if the argmax multiplier is
in the range; the question is what happens beyond.
Run: python 04-computation/experiments/collatz_five_mirrors_multiplier_families_20260929.py [HMAX] [UMAX]   (~10 min)
"""
from __future__ import annotations

import cmath
import math
import sys
from fractions import Fraction


def family_max(u: int, HMAX: int, AMAX: int = 40, JMAX_EXTRA: int = 20):
    JMAX = 2 * HMAX + JMAX_EXTRA
    weights = [2.0 ** (-a) for a in range(1, AMAX + 1)]
    prev = {j: 1.0 + 0j for j in range(-AMAX * HMAX - AMAX, JMAX + 1)}
    out = []
    for n in range(1, HMAX + 1):
        mod = 3 ** n
        lo = -AMAX * (HMAX - n)
        pw = {j: (u * pow(2, j, mod)) % mod for j in range(lo - AMAX, JMAX + 1)}
        cur = {}
        for j in range(lo, JMAX + 1):
            s = 0j
            for a in range(1, AMAX + 1):
                r = pw[j - a]
                s += weights[a - 1] * cmath.exp(2j * math.pi * float(Fraction(r, mod))) * prev.get(j - a, 0j)
            cur[j] = s
        best = max(cur, key=lambda j: abs(cur[j]))
        out.append((n, abs(cur[best]), best))
        prev = cur
    return out


if __name__ == "__main__":
    HMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    UMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 49
    us = [u for u in range(1, UMAX + 1, 2) if u % 3]
    res = {u: family_max(u, HMAX) for u in us}
    print(f"multipliers u = {us}; levels h <= {HMAX}")
    print("h: max over all families, its multiplier u* (and exponent offset s-h), the pure-family value, ratio pure/max")
    for i in range(HMAX):
        h = i + 1
        best_u = max(us, key=lambda u: res[u][i][1])
        mx = res[best_u][i][1]
        pure = res[1][i][1]
        print(f"   h={h:3d}: max {mx:.6e} at u={best_u} (s-h={res[best_u][i][2]-h:+d}); pure family {pure:.6e}; pure/max = {pure/mx:.4f}"
              + ("" if best_u == 1 else "   <-- pure family not maximal"))
    print("DONE")
