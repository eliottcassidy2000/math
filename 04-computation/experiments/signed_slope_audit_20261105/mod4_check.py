#!/usr/bin/env python3
"""Independent numeric check of THM-4882 (2) Lemma 1 and (3)(iii) for many even-denominator slopes, n < 2^NBITS.

For each coprime (a, b), a odd, b even, 1/2 < a/b < 1:
  4a <= 3b: every positive value of F(n) = a n - b g(n) (n odd >= 3) must be = -(a+b) mod 4;
            also report the classes of negative values, and the odd-m branch sign (must be negative).
  4a > 3b : every value below -N_0 (N_0 = a + ab/(4a-3b)) must be = 2 - a (4 | b) or -a (b = 2 mod 4) mod 4;
            report the largest even-m value that is negative (to see how sharp N_0 is) and the classes above -N_0.
Also: the "misses infinitely many" conclusion is checked as: which odd targets in [-2000, 2000] are missed at this n range.
"""
import sys
import numpy as np
from math import gcd
from fractions import Fraction

NBITS = int(sys.argv[1]) if len(sys.argv) > 1 else 24
top = 1 << NBITS
g = np.zeros(top + 1, dtype=np.int64); g[1] = 1
# fill in increasing order: g(2n) = g(n), g(2n-1) = g(n) + n - (n mod 2)
j = 1
while (1 << j) <= top:
    lo = 1 << j; hi = min(1 << (j + 1), top + 1)
    ev = np.arange(lo, hi, 2, dtype=np.int64); g[ev] = g[ev // 2]
    od = np.arange(lo + 1, hi, 2, dtype=np.int64); mm = (od + 1) // 2; g[od] = g[mm] + mm - (mm % 2)
    j += 1
# spot-check against the plain recurrence
def gr(n):
    if n == 1: return 1
    if n % 2 == 0: return gr(n // 2)
    m = (n + 1) // 2; return gr(m) + m - (m % 2)
import random
rs = random.Random(1)
assert all(g[n] == gr(n) for n in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,31,32,33] + [rs.randrange(1, top + 1) for _ in range(2000)])
nodd = np.arange(3, top, 2, dtype=np.int64); godd = g[3:top:2]
# Lemma 1
assert np.all(g[1:] % 2 == 1), "g not odd"
assert np.all(g[1:] <= np.arange(1, top + 1)), "g > n"
assert np.all(2 * godd >= nodd + 3), "lower bound (n+3)/2 fails"
# tightness of the lower bound
tight = nodd[2 * godd == nodd + 3]
print(f"Lemma 1 verified for n <= 2^{NBITS}: g odd, g <= n, g(n) >= (n+3)/2 for odd n >= 3; equality at n = {tight[:8].tolist()} ... ({len(tight)} cases: all of the form 2^j - 1? {all((t + 1) & t == 0 for t in tight.tolist())})")

slopes = []
for b in range(4, 26, 2):
    for a in range(3, b, 2):
        if gcd(a, b) == 1 and 2 * a > b: slopes.append((a, b))
for extra in [(17, 18), (19, 20), (11, 20), (13, 20), (17, 24), (15, 16), (9, 10), (5, 8), (3, 4), (23, 24), (19, 24), (13, 24), (25, 26), (15, 26), (31, 32), (17, 32)]:
    if extra not in slopes: slopes.append(extra)
viol = 0
for (a, b) in slopes:
    F = a * nodd - b * godd
    pos = F[F > 0]; neg = F[F < 0]
    m = (nodd + 1) // 2
    oddm = F[m % 2 == 1]; evenm = F[m % 2 == 0]
    if 4 * a <= 3 * b:
        cls = set((pos % 4).tolist()); pred = (-(a + b)) % 4
        okc = cls <= {pred}
        oddm_max = oddm.max()
        line = f"{a}/{b} (<= 3/4): positive classes mod 4 {sorted(cls)} predicted {{{pred}}} {'OK' if okc else 'VIOLATION'}; odd-m branch max {oddm_max} (must be < 0: {'OK' if oddm_max < 0 else 'VIOLATION'}; bound -a-b/2 = {-a - b//2}); negative classes {sorted(set((neg % 4).tolist()))}"
        if not okc or oddm_max >= 0: viol += 1
    else:
        N0 = Fraction(a) + Fraction(a * b, 4 * a - 3 * b)
        below = F[F < -N0]
        cls = set((below % 4).tolist()); pred = (2 - a) % 4 if b % 4 == 0 else (-a) % 4
        okc = cls <= {pred}
        evenm_neg = evenm[evenm < 0]
        evenm_min = evenm_neg.min() if len(evenm_neg) else None
        above_cls = sorted(set((neg[neg >= -N0] % 4).tolist()))
        line = f"{a}/{b} (> 3/4): N_0 = {N0} = {float(N0):.2f}; classes below -N_0 {sorted(cls)} predicted {{{pred}}} {'OK' if okc else 'VIOLATION'}; most negative even-m value {evenm_min} (must be > -N_0: {'OK' if evenm_min is None or evenm_min > -N0 else 'VIOLATION'}); negative classes in [-N_0, 0) {above_cls}; positive classes {sorted(set((pos % 4).tolist()))}"
        if not okc or (evenm_min is not None and evenm_min <= -N0): viol += 1
    # missed odd targets in [-2000, 2000] (small-witness)
    att = set(F[(F >= -2000) & (F <= 2000)].tolist())
    missed = [t for t in range(-1999, 2000, 2) if t not in att]
    line += f"; odd targets in [-2000,2000] missed with n < 2^{NBITS}: {len(missed)} of 2000"
    print(line)
print(f"slopes checked: {len(slopes)}; violations: {viol}")
