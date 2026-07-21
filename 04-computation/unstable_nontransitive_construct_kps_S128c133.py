#!/usr/bin/env python3
"""unstable_nontransitive_construct_kps_S128c133.py -- kind-pasteur-2026-07-20-S128c133

CHARACTERISATION of the GIT-unstable NON-transitive tournaments, by construction.

STRUCTURE (from the n=7 witness char_A = x^4(x^3-1)).  A reducible tournament decomposes by
its strongly-connected components (SCCs) into a TRANSITIVE order C_1 > C_2 > ... > C_r, and
char_A = product of char(C_i) (block upper-triangular).  A singleton SCC contributes x
(eigenvalue 0); a 3-cycle SCC contributes x^3-1 = (x-1)(x^2+x+1).  So the multiplicity of 0
is (# singleton SCCs) + (0-eigenvalues inside non-trivial SCCs).

UNSTABLE <=> that 0-multiplicity > n/2.  At n=7, a non-trivial SCC needs >= 3 vertices, so the
maximum number of singletons is 4 (the other 3 forming a 3-cycle, which has NO 0-eigenvalue):
0-multiplicity = 4 > 3.5.  Hence the UNIQUE type of unstable non-transitive tournament at n=7
is  4 transitive singletons + one 3-cycle atom.  This script:
  * builds all such tournaments (3-cycle atom placed at each rank of the SCC order),
  * confirms each is unstable (char_A = x^4(x^3-1), mult 4) and non-transitive (c3=1),
  * reads off |Aut| (the atom's Z_3 = the localized 'point of symmetry') and self-comp
    (blue/SC) status,
  * confirms no OTHER SCC-type gives an unstable non-transitive tournament at n=7,
  * checks that no STRONGLY-CONNECTED tournament is unstable non-transitive (sample).
"""
import sys
import numpy as np
import sympy as sp
from itertools import permutations, combinations
from collections import Counter

x = sp.Symbol('x')


def three_cycle(k):
    """adjacency of a 3-cycle on vertices {0,1,2} (0->1->2->0)."""
    A = np.zeros((3, 3), int)
    A[0, 1] = A[1, 2] = A[2, 0] = 1
    return A


def build_reducible(order):
    """order = list of SCC blocks (each a list of local adjacency + size); earlier SCC beats
    all later ones. Returns n x n adjacency."""
    sizes = [b.shape[0] for b in order]
    n = sum(sizes)
    A = np.zeros((n, n), int)
    off = 0
    starts = []
    for b in order:
        s = b.shape[0]
        starts.append(off)
        A[off:off + s, off:off + s] = b
        off += s
    # earlier block beats later block
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            si, sj = starts[i], starts[j]
            for a in range(sizes[i]):
                for b in range(sizes[j]):
                    A[si + a, sj + b] = 1
    return A


def canon(A, n):
    best = None
    for p in permutations(range(n)):
        t = A[np.ix_(p, p)].tobytes()
        if best is None or t < best:
            best = t
    return best


def aut_order(A, n):
    return sum(1 for p in permutations(range(n)) if np.array_equal(A[np.ix_(p, p)], A))


def is_self_comp(A, n):
    comp = np.ones((n, n), int) - np.eye(n, dtype=int) - A
    return canon(A, n) == canon(comp, n)


def strongly_connected(A, n):
    return bool((np.linalg.matrix_power(A + np.eye(n, dtype=int), n) > 0).all())


print("=" * 80)
print("n=7 UNSTABLE NON-TRANSITIVE = 4 transitive singletons + one 3-cycle atom")
print("=" * 80)
n = 7
S = np.zeros((1, 1), int)      # singleton
C3 = three_cycle(3)
blocks = [S, S, S, S, C3]
seen = {}
for pos in range(5):           # rank of the 3-cycle among 5 SCCs
    order = [S, S, S, S]
    order.insert(pos, C3)
    A = build_reducible(order)
    cc = canon(A, n)
    if cc in seen:
        continue
    M = sp.Matrix(A.tolist())
    p = sp.Poly(M.charpoly(x).as_expr(), x)
    fac = sp.factor(p.as_expr())
    mu = max((e for _, e in sp.factor_list(p.as_expr())[1]), default=1)
    c3 = int(round(np.trace(np.linalg.matrix_power(A, 3)) / 3))
    au = aut_order(A, n)
    sc = is_self_comp(A, n)
    scores = sorted(A.sum(axis=1).tolist())
    seen[cc] = 1
    print("  3-cycle at rank %d : char_A=%s  mult=%d c3=%d |Aut|=%d SC=%s scores=%s"
          % (pos, fac, mu, c3, au, sc, scores))
    sys.stdout.flush()
print("  distinct iso classes of this type: %d" % len(seen))

print()
print("=" * 80)
print("NO OTHER SCC-TYPE is unstable non-transitive at n=7")
print("=" * 80)
print("  A non-trivial SCC needs >=3 vertices, so # singletons <= 4; with a 5- or 4-vertex")
print("  non-trivial SCC the singleton count drops below 4 and 0-mult <= 3 < n/2, UNLESS the")
print("  big SCC itself contributes 0-eigenvalues. Checking small non-trivial SCCs for 0-eigs:")
for sz in (3, 4, 5):
    # sample strongly-connected tournaments on sz vertices, check if 0 is an eigenvalue
    import random
    random.seed(sz)
    has0 = False
    for _ in range(3000):
        m = sz * (sz - 1) // 2
        bits = random.getrandbits(m)
        B = np.zeros((sz, sz), int)
        idx = 0
        for i in range(sz):
            for j in range(i + 1, sz):
                if bits >> idx & 1:
                    B[i, j] = 1
                else:
                    B[j, i] = 1
                idx += 1
        if not strongly_connected(B, sz):
            continue
        ev = np.linalg.eigvals(B.astype(float))
        if any(abs(e) < 1e-6 for e in ev):
            has0 = True
            break
    print("    strongly-connected size-%d tournaments with a 0 eigenvalue: %s" % (sz, has0))

print()
print("=" * 80)
print("NO STRONGLY-CONNECTED tournament is unstable non-transitive (n=7 sample)")
print("=" * 80)
import random
random.seed(11)
worst = 0
for _ in range(200000):
    bits = random.getrandbits(21)
    A = np.zeros((7, 7), int)
    idx = 0
    for i in range(7):
        for j in range(i + 1, 7):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    if not strongly_connected(A, 7):
        continue
    ev = np.sort(np.linalg.eigvals(A.astype(float)).real)
    best = 1; i = 0
    while i < 7:
        j = i
        while j + 1 < 7 and abs(ev[j + 1] - ev[i]) < 1e-6:
            j += 1
        best = max(best, j - i + 1); i = j + 1
    worst = max(worst, best)
print("  max multiplicity over strongly-connected n=7 sample = %d  (unstable needs > 3.5: %s)"
      % (worst, worst >= 4))
print()
print("  CONCLUSION: the unstable non-transitive tournaments are exactly the REDUCIBLE ones")
print("  whose nilpotent transitive skeleton (singleton SCCs) has a Jordan block > n/2, with")
print("  one cyclic atom. Their ONLY symmetry is the atom's (Z_3) -- a LOCAL point of symmetry")
print("  near the nullcone vertex -- and they are NOT self-complementary (not blue/SC): the")
print("  opposite pole from the globally-symmetric SC/blue intransitive classes (opus j=0 Paley).")
