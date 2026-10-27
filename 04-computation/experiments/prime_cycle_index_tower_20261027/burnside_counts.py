"""Tournament classes A(n) and self-converse classes SC(n) by Burnside over cycle types (exact).

For a permutation s with cycle type lam, the s-orbits on unordered pairs each have a length L and a
flag 'swap' (s^L maps {u,v} to itself exchanging u and v).
  s(T) = T     : consistent iff no orbit swaps;               count 2^(#orbits)
  s(T) = T^op  : consistent iff every orbit has (L even, no swap) or (L odd, swap)
A(n)  = sum_lam [consistent] 2^orbits / z_lam   (Davis 1954)
SC(n) = sum_lam [anti-consistent] 2^orbits / z_lam   (classes T with T ~ T^op: 2*#orbits(S_n x Z_2) - #orbits(S_n)).
"""
from fractions import Fraction
from math import factorial, gcd
from collections import Counter
import sys

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield []; return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k): yield [k] + rest

def z(lam):
    c = Counter(lam); r = 1
    for l, m in c.items(): r *= l ** m * factorial(m)
    return r

def pair_orbits(lam):
    """list of (L, swap) for the orbits of a permutation of cycle type lam on unordered pairs"""
    res = []
    # within a cycle of length l: pairs at distance d, 1 <= d <= l//2
    for l in lam:
        for d in range(1, l // 2 + 1):
            if 2 * d == l: res.append((l // 2, True))      # antipodal pairs: orbit of size l/2, s^(l/2) swaps
            else: res.append((l, False))
    # between cycles of lengths a, b: gcd(a,b) orbits of length lcm(a,b), never swapping
    for i in range(len(lam)):
        for j in range(i + 1, len(lam)):
            a, b = lam[i], lam[j]; g = gcd(a, b)
            for _ in range(g): res.append((a * b // g, False))
    return res

def counts(n):
    A = Fraction(0); SC = Fraction(0)
    for lam in partitions(n):
        orb = pair_orbits(lam)
        w = Fraction(1, z(lam)) * 2 ** len(orb)
        if all(not s for (L, s) in orb): A += w
        if all((L % 2 == 0 and not s) or (L % 2 == 1 and s) for (L, s) in orb): SC += w
    assert A.denominator == 1 and SC.denominator == 1
    return int(A), int(SC)

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    for n in range(1, NMAX + 1):
        a, sc = counts(n)
        print(n, a, sc)
