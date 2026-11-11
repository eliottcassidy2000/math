#!/usr/bin/env python3
"""Arithmetic of the multiplicative monoid M = {odd positive integers} minus {7, 21} (the conjectured set of
Hamiltonian-path counts).  Atoms (KNOWN: THM-4094, WEIGHTED_ATLAS (9)): the odd primes other than 7, the numbers 7p
with p an odd prime other than 3 (49 included), 63 and 343.  This script lists the atoms below a bound by brute
force, computes the shortest and longest factorisations of every element into atoms, and reports the elasticity
(largest ratio longest/shortest).  Claim of the note: the elasticity is 5/2, attained at 63^2 = 3^4 * 49.
Usage: hmonoid_elasticity.py [bound, default 300000]"""
import sys
from fractions import Fraction as Fr
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
inM = [False] * (N + 1)
for m in range(1, N + 1, 2): inM[m] = True
inM[7] = inM[21] = False
lo = [0] * (N + 1); hi = [0] * (N + 1); atoms = []
for m in range(3, N + 1, 2):
    if not inM[m]: continue
    best_lo = None; best_hi = None
    d = 3
    while d * d <= m:
        if m % d == 0:
            for a in (d, m // d):
                b = m // a
                if inM[a] and inM[b] and a > 1 and b > 1:
                    l = lo[a] + lo[b]; h = hi[a] + hi[b]
                    best_lo = l if best_lo is None else min(best_lo, l)
                    best_hi = h if best_hi is None else max(best_hi, h)
        d += 2
    if best_lo is None: atoms.append(m); lo[m] = hi[m] = 1
    else: lo[m] = best_lo; hi[m] = best_hi
def isprime(k):
    if k < 2: return False
    d = 2
    while d * d <= k:
        if k % d == 0: return False
        d += 1
    return True
pred = sorted([p for p in range(3, N + 1, 2) if isprime(p) and p != 7] + [7 * p for p in range(5, N // 7 + 1, 2) if isprime(p)] + [x for x in (63, 343) if x <= N])
print("bound %d: %d atoms; equal to the predicted list (primes != 7, 7p with p != 3, 63, 343): %s" % (N, len(atoms), atoms == pred))
print("composite atoms below 400:", [a for a in atoms if not isprime(a) and a < 400])
best = max((Fr(hi[m], lo[m]), -m) for m in range(3, N + 1, 2) if inM[m])
print("largest ratio longest/shortest factorisation below the bound: %s, first attained at %d (shortest %d, longest %d)" % (best[0], -best[1], lo[-best[1]], hi[-best[1]]))
print("elements with ratio 5/2:", [m for m in range(3, N + 1, 2) if inM[m] and Fr(hi[m], lo[m]) == Fr(5, 2)][:6])
print("smallest element with two factorisations of different lengths:", next(m for m in range(3, N + 1, 2) if inM[m] and hi[m] != lo[m]))
