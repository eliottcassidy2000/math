# Exact E(L,k) = Lyn(L,k)/|2^L-3^k| for 12 <= L <= LMAX, all k >= L/2; report the top shapes.
import sys
from math import comb, gcd
from fractions import Fraction
from sympy import mobius, divisors
LMAX = int(sys.argv[1])
def lyn(L, k):
    g = gcd(L, k) if k else L
    return sum(int(mobius(e)) * comb(L // e, k // e) for e in divisors(g)) // L
best = []
for L in range(12, LMAX + 1):
    twoL = 1 << L
    # only k near L*log_3 2 can matter; scan a window, but use exact ints
    k0 = int(L / 1.5849625007211562)
    for k in range(max((L + 1) // 2, k0 - 3), min(L, k0 + 4) + 1):
        m = abs(twoL - 3**k)
        # compare via floats of logs to avoid huge Fractions
        from math import log
        lE = (lyn(L, k).bit_length() - m.bit_length())
        if lE > -12:   # E > ~2^-12
            E = lyn(L, k) / m
            best.append((E, L, k))
best.sort(reverse=True)
print("top shapes with L>=12, L<=%d:" % LMAX)
for E, L, k in best[:12]:
    print("  (%d,%d) E=%.5f" % (L, k, E))
