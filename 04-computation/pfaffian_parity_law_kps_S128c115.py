#!/usr/bin/env python3
"""pfaffian_parity_law_kps_S128c115.py -- kind-pasteur-2026-07-20-S128c115

WHY THE SKEW CHAR-POLY COEFFICIENTS ARE CONGRUENT -- and it is the odd/even mechanism.

The n = 7 search found only 11 distinct skew characteristic polynomials, with attained
coefficients

    x^3 : 35, 67, 83, 99, 115, 131, 147      all  == 3  (mod 16)
    x^1 : 7, 23, 39, 55, 71, 119, 135, ...   all  == 7  (mod 16)

That is not an accident, and the reason is exactly the odd/even theme:

    for a real skew-symmetric S, det S = Pf(S)^2 in EVEN order, and det S = 0 in ODD
    order.  So in the expansion of char(S) only EVEN principal minors survive, and each
    equals a squared Pfaffian:

        coefficient of x^{n-2k} in char(S)  =  sum over 2k-subsets A of Pf(S_A)^2 .

    Odd-order minors contribute NOTHING.  That single fact is what forces char(S) to be
    parity-pure (odd for n odd, even for n even) -- the spectral form of "odd functions
    have no even part".

CONSEQUENCES, all checked here for tournaments (S_ij = +-1 off the diagonal):

  k = 1: Pf of a 2x2 is S_12, so the x^{n-2} coefficient is sum over PAIRS of 1
         = C(n,2).  Carries no information -- which is why 21 appeared in both the
         owner's polynomial and Paley's.
  k = 2: Pf of a 4x4 is S12*S34 - S13*S24 + S14*S23, a sum of three +-1 terms, so
         Pf is +-1 or +-3 and Pf^2 in {1,9}.  Hence
             c_{n-4} = C(n,4) + 8*k4,   k4 = #{4-subsets with |Pf| = 3}.
         The observed congruence mod 16 then says k4 is EVEN.  Tested.
  k = 3: Pf of a 6x6 is a sum of 15 +-1 terms, hence ODD, so Pf^2 == 1 (mod 8) and
             c_{n-6} == C(n,6)  (mod 8).
         The observed mod 16 needs more than parity; tested rather than assumed.
"""
import sys
import random
from itertools import combinations
import numpy as np

random.seed(29)
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 4000


def from_bits(bits, n):
    adj = [0] * n
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            idx += 1
    return adj


def skew(adj, n):
    S = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if i != j:
                S[i, j] = 1 if (adj[i] >> j) & 1 else -1
    return S


def charpoly(A, n):
    I = np.eye(n, dtype=np.int64)
    M = I.copy()
    cs = [1]
    for k in range(1, n + 1):
        AM = A.dot(M)
        c = -int(np.trace(AM)) // k
        cs.append(c)
        M = AM + c * I
    return cs


def pf(S, idx):
    """Pfaffian of the principal submatrix on idx (even size), by recursion."""
    m = len(idx)
    if m == 0:
        return 1
    if m == 2:
        return int(S[idx[0], idx[1]])
    tot = 0
    i0 = idx[0]
    for t in range(1, m):
        j = idx[t]
        rest = [idx[u] for u in range(1, m) if u != t]
        tot += ((-1) ** (t - 1)) * int(S[i0, j]) * pf(S, rest)
    return tot


print("=" * 78)
print("MINOR EXPANSION: coeff of x^{n-2k} in char(S)  ==  sum over 2k-subsets Pf^2")
print("=" * 78)
for n in (5, 6, 7):
    m = n * (n - 1) // 2
    bad = 0
    for _ in range(40):
        adj = from_bits(random.getrandbits(m), n)
        S = skew(adj, n)
        cs = charpoly(S, n)
        for k in range(1, n // 2 + 1):
            lhs = cs[2 * k]
            rhs = sum(pf(S, list(A)) ** 2 for A in combinations(range(n), 2 * k))
            if lhs != rhs:
                bad += 1
    print("  n = %d : mismatches over 40 tournaments and all k : %d  -> %s"
          % (n, bad, "IDENTITY HOLDS" if bad == 0 else "IDENTITY FAILS"))

print()
print("=" * 78)
print("THE 4-SUBSET LAW:  c_{n-4} = C(n,4) + 8*k4,  and is k4 EVEN?")
print("=" * 78)
for n in (6, 7, 8):
    m = n * (n - 1) // 2
    okform = True
    par = {0: 0, 1: 0}
    vals = set()
    for _ in range(NS if n < 8 else NS // 4):
        adj = from_bits(random.getrandbits(m), n)
        S = skew(adj, n)
        cs = charpoly(S, n)
        k4 = sum(1 for A in combinations(range(n), 4) if abs(pf(S, list(A))) == 3)
        from math import comb
        if cs[4] != comb(n, 4) + 8 * k4:
            okform = False
        par[k4 % 2] += 1
        vals.add(cs[4] % 16)
    from math import comb
    print("  n = %d : c_{n-4} = C(n,4) + 8*k4 always : %-5s   C(n,4) = %d"
          % (n, okform, comb(n, 4)))
    print("          k4 parity: even %d, odd %d   -> k4 always even : %s"
          % (par[0], par[1], par[1] == 0))
    print("          c_{n-4} mod 16 attained : %s" % sorted(vals))

print()
print("=" * 78)
print("THE 6-SUBSET LAW:  Pf of a 6x6 +-1 skew matrix is ODD, so Pf^2 == 1 (mod 8)")
print("=" * 78)
for n in (6, 7):
    m = n * (n - 1) // 2
    allodd = True
    vals = set()
    from math import comb
    for _ in range(600):
        adj = from_bits(random.getrandbits(m), n)
        S = skew(adj, n)
        for A in combinations(range(n), 6):
            if pf(S, list(A)) % 2 == 0:
                allodd = False
        cs = charpoly(S, n)
        vals.add((cs[6] - comb(n, 6)) % 8)
    print("  n = %d : every 6x6 Pfaffian odd : %-5s ; (c_{n-6} - C(n,6)) mod 8 : %s"
          % (n, allodd, sorted(vals)))

print()
print("=" * 78)
print("READING")
print("=" * 78)
print("  The parity of char(S) is not a curiosity about tournaments -- it is the")
print("  Pfaffian.  Odd-order skew minors vanish identically, so only even-order ones")
print("  enter, and each enters as a SQUARE.  Squares of integers are what produce the")
print("  congruences: |Pf| in {1,3} at size 4 gives steps of 8, and oddness of Pf at")
print("  size 6 gives 1 mod 8.  'Odd function' and 'even subgraph' meet here: the")
print("  even-size principal minors are the only ones that survive, and they survive")
print("  squared.")
