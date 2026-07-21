#!/usr/bin/env python3
"""unstable_nontransitive_structure_kps_S128c133.py -- kind-pasteur-2026-07-20-S128c133

CHARACTERISE the GIT-unstable NON-transitive tournaments (THM-1825): those whose char_A has a
root of multiplicity > n/2 (a Jordan block, Lemma A/B) but that are not transitive.  They first
appear at n=7.  Enumerate them exhaustively at n=7, dedup by iso class, and read off:
  * reducibility (strongly connected?  a dominating / dominated vertex?),
  * automorphism group order (a "point of symmetry"?),
  * self-complementary (= blue / grid-symmetric / SC)?,
  * the exact char_A factorisation (the partial-collision j-stratum near the nullcone vertex).

Frame: the unstable non-transitive sit NEAR the nullcone vertex x^n (transitive) -- partially
collided roots -- opposite to the symmetric/blue (SC, j=0 equianharmonic) intransitive pole
where opus-S434 put the MOST intransitive (Paley).  Test which pole they live at.
"""
import sys
import numpy as np
import sympy as sp
from itertools import permutations
from collections import Counter

x = sp.Symbol('x')


def from_bits(bits, n):
    A = np.zeros((n, n), dtype=np.int64)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    return A


def numpy_maxmult(A, n):
    ev = np.sort(np.linalg.eigvals(A.astype(float)).real)
    best = 1; i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(ev[j + 1] - ev[i]) < 1e-6:
            j += 1
        best = max(best, j - i + 1); i = j + 1
    return best


def canon(A, n):
    best = None
    for p in permutations(range(n)):
        B = A[np.ix_(p, p)]
        t = B.tobytes()
        if best is None or t < best:
            best = t
    return best


def aut_order(A, n):
    cnt = 0
    for p in permutations(range(n)):
        if np.array_equal(A[np.ix_(p, p)], A):
            cnt += 1
    return cnt


def is_self_comp(A, n):
    """self-complementary: A^op = J-I-A is isomorphic to A (blue/SC)."""
    comp = (np.ones((n, n), int) - np.eye(n, dtype=int) - A)
    ca = canon(A, n)
    cc = canon(comp, n)
    return ca == cc


def strongly_connected(A, n):
    reach = (np.linalg.matrix_power(A + np.eye(n, dtype=np.int64), n) > 0)
    return bool(reach.all())


n = 7
m = n * (n - 1) // 2
print("=" * 82)
print("EXHAUSTIVE n=7 scan for unstable non-transitive tournaments (2^%d)" % m)
print("=" * 82)
sys.stdout.flush()
seen = set()
classes = []            # (canon, A)
flagged = 0
for bits in range(1 << m):
    A = from_bits(bits, n)
    if numpy_maxmult(A, n) < 4:          # need mult > n/2 = 3.5, i.e. >= 4
        continue
    flagged += 1
    # exact check
    M = sp.Matrix(A.tolist())
    p = sp.Poly(M.charpoly(x).as_expr(), x)
    mu = max((e for _, e in sp.factor_list(p.as_expr())[1]), default=1)
    if mu < 4:                            # numpy clustering artifact
        continue
    c3 = int(round(np.trace(np.linalg.matrix_power(A, 3)) / 3))
    if c3 == 0:                           # transitive -- excluded
        continue
    cc = canon(A, n)
    if cc in seen:
        continue
    seen.add(cc)
    classes.append((A.copy(), p, c3))
print("  numpy-flagged: %d ; genuine unstable non-transitive ISO CLASSES: %d"
      % (flagged, len(classes)))
sys.stdout.flush()

print()
print("=" * 82)
print("STRUCTURE of each unstable non-transitive iso class")
print("=" * 82)
print("  %-28s %-6s %-5s %-6s %-5s %-6s %s"
      % ("char_A (factored)", "c3", "scc?", "|Aut|", "SC?", "scores", "eig-mult"))
red_all = True
sc_any = False
aut_orders = Counter()
for A, p, c3 in classes:
    fac = sp.factor(p.as_expr())
    scc = strongly_connected(A, n)
    au = aut_order(A, n)
    sc = is_self_comp(A, n)
    scores = sorted(A.sum(axis=1).tolist())
    mu = max((e for _, e in sp.factor_list(p.as_expr())[1]), default=1)
    red_all = red_all and (not scc)
    sc_any = sc_any or sc
    aut_orders[au] += 1
    print("  %-28s %-6d %-5s %-6d %-5s %-14s mult=%d"
          % (str(fac)[:27], c3, scc, au, sc, str(scores), mu))
    sys.stdout.flush()

print()
print("=" * 82)
print("VERDICT")
print("=" * 82)
print("  all REDUCIBLE (not strongly connected): %s" % red_all)
print("  any self-complementary (blue/SC): %s" % sc_any)
print("  |Aut| distribution: %s" % dict(aut_orders))
print()
print("  Reading: unstable non-transitive tournaments sit NEAR the nullcone vertex (a root")
print("  of multiplicity > n/2 = a big transitive/nilpotent Jordan block), and this test")
print("  shows whether they are reducible (a dominating/dominated vertex + a small cyclic")
print("  atom) and whether they carry the 'point of symmetry' (nontrivial Aut) -- vs the")
print("  blue/SC symmetric-intransitive pole where opus-S434 put the MOST intransitive.")
