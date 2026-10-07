"""Rational cycles of the Terras map with shape (L,k) = Ellison exceptions and (11,7):
x_w = c_w/(2^L - 3^k) over all Lyndon words w with k ones; tally reduced denominators
(a reduced denominator d means an integer cycle of the 3x+d map of the same shape, Lagarias 1990)."""
from math import gcd
from collections import Counter
from sympy import factorint, divisors
def lyndon_words_k(L, k):
    # FKM with fixed number of ones (generate all Lyndon words, filter by count) -- fine for L<=27
    a = [0] * (L + 1)
    out = []
    def gen(t, p, ones):
        if ones > k or ones + (L - t + 1) < k: return
        if t > L:
            if p == L and ones == k: out.append(a[1:L+1][:])
            return
        a[t] = a[t - p]; gen(t + 1, p, ones + a[t])
        if a[t - p] == 0:
            a[t] = 1; gen(t + 1, t, ones + 1)
    gen(1, 1, 0)
    return out
def cw(w):
    d = 0
    for t, b in enumerate(w):
        if b: d = 3 * d + (1 << t)
    return d
for (L, k) in [(11, 7), (13, 8), (14, 9), (16, 10), (19, 12), (27, 17)]:
    m = 2**L - 3**k
    W = lyndon_words_k(L, k)
    den = Counter(); bestfrac = None
    for w in W:
        c = cw(w); g = gcd(c, abs(m)); den[abs(m) // g] += 1
        # distance of x_w to nearest integer, in units of 1/|m|
        r = c % abs(m); dist = min(r, abs(m) - r)
        if bestfrac is None or dist < bestfrac[0]: bestfrac = (dist, c, w)
    small = sorted((d, n) for d, n in den.items() if d < abs(m))
    print("(L,k)=(%d,%d) m=%d=%s  Lyndon words=%d  denominators<|m| realised: %s" % (L, k, m, dict(factorint(abs(m))), len(W), small))
    dist, c, w = bestfrac
    print("     closest-to-integer x_w: c_w=%d, x_w=%.6f (dist %d/|m|), word=%s" % (c, c / m, dist, "".join(map(str, w))))
