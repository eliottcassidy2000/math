#!/usr/bin/env python3
"""Tables of the parity races by cycle length from the outputs of cycrace.c (out/cyc_n.txt, n = 3..10).
Each line: K signed weight, K a bitmask over the odd lengths 3, 5, 7, 9; the value is E (-1)^(sum_{k in K} c_k) over
labelled tournaments.  Also prints the same characters in the subset basis
    Y_k = sum over k-subsets S of c_odd(T[S])  =  (p_k(T) - C(n,k))/2  (mod 2),   p_k = number of k-vertex paths,
using  Y_k = sum_l C(n - l, k - l) c_l  (mod 2), and checks the identities recorded in the note."""
from fractions import Fraction
from math import comb
import os
here = os.path.dirname(os.path.abspath(__file__))
L = [3, 5, 7, 9]
def load(n):
    d = {}
    for line in open(os.path.join(here, 'out', 'cyc_%d.txt' % n)):
        K, s, w = map(int, line.split()); d[K] = Fraction(s, w)
    return d
def name(K, nl): return "+".join("c%d" % L[b] for b in range(nl) if K >> b & 1) or "1"
E = {}
for n in range(3, 11):
    d = load(n); nl = len([k for k in range(3, n + 1, 2)]); E[n] = (d, nl)
    print("n = %2d: " % n + " | ".join("%s: %s" % (name(K, nl), d[K]) for K in sorted(d) if K))
print()
print("Checks")
for n in range(3, 11):
    d, nl = E[n]; assert d[1] == Fraction(1, 2 ** ((n - 1) // 2)), n
print("  triangle law E(-1)^c3 = 2^(-floor((n-1)/2)) for n = 3..10: OK")
full = {n: E[n][0][(1 << E[n][1]) - 1] for n in E}
print("  full race D_n, n = 3..10:", ", ".join(str(full[n]) for n in sorted(full)))
print()
print("Subset basis: Y_k = sum_l C(n-l, k-l) c_l mod 2 for odd k")
for n in range(5, 11):
    d, nl = E[n]
    rows = []
    for k in [x for x in L[:nl]]:
        mask = 0
        for b in range(nl):
            if L[b] <= k and comb(n - L[b], k - L[b]) % 2: mask |= 1 << b
        rows.append((k, mask))
    print("  n = %2d: " % n + " ; ".join("Y_%d = %s -> E = %s" % (k, name(m, nl), d[m]) for k, m in rows))
print()
print("Recorded coincidences (FINITE-EXACT) and failed extrapolations")
d6, d8, d10, d9 = E[6][0], E[8][0], E[10][0], E[9][0]
print("  n = 6:  E[c3+c5] = E[c3]:", d6[3] == d6[1], ";  E[c5] = E[c3]^2:", d6[2] == d6[1] ** 2)
print("  n = 8:  E[c3+c5+c7] = E[c3+c7]:", d8[7] == d8[5])
print("  n = 10: E[c3+c5+c7+c9] = E[c3+c7+c9] (R5):", d10[15] == d10[13], ";  E[c3+c5] = E[c3]:", d10[3] == d10[1], ";  E[c5] = E[c3]^2:", d10[2] == d10[1] ** 2)
print("  n = 10: E[c7] = E[c3+c5] E[c3+c5+c7]:", d10[4] == d10[3] * d10[7], ";  E[c3+c7] = E[c3+c5] E[c5+c7]:", d10[5] == d10[3] * d10[6])
print("  pair law n = 7 -> 8 for c3, c5, c3+c5:", [E[7][0][K] == d8[K] for K in (1, 2, 3)])
print("  pair law n = 9 -> 10 for c3, c5, c3+c5, c7, c3+c7, c5+c7, c3+c5+c7:", [d9[K] == d10[K] for K in (1, 2, 3, 4, 5, 6, 7)])
