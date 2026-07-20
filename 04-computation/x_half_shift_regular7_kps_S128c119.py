#!/usr/bin/env python3
"""x_half_shift_regular7_kps_S128c119.py -- kind-pasteur-2026-07-20-S128c119

THE OWNER'S x(x^2+7)(x^4+14x^2+17), FOUND.

THM-1555 says the half-dictionary is A = (J - I + S)/2, i.e. S = 2A - J + I.  On the
all-ones complement of a REGULAR tournament J vanishes, so the eigenvalues correspond by

    x = 2*lambda + 1

-- literally the inverse of the owner's map x |-> (1+x)/2, undoubled.  So substituting
x = 2 lambda + 1 into the characteristic polynomial of A should turn the repo's
adjacency spectra into the SKEW spectra, and the skew spectra are purely imaginary,
which means the result must be an EVEN polynomial (times the Perron factor x).

That is exactly the owner's even/odd axis, and it predicts the shape

    char_S(x) = x * (even polynomial in x).

This enumerates every regular tournament on 7 vertices up to isomorphism and prints
char_S in the x coordinate, to see which one is x(x^2+7)(x^4+14x^2+17).
"""
from itertools import combinations, permutations
import numpy as np
import sympy as sp

x = sp.Symbol('x')
n = 7
pairs = list(combinations(range(n), 2))

def canon(A):
    best = None
    for p in permutations(range(n)):
        B = A[np.ix_(p, p)]
        t = B.tobytes()
        if best is None or t < best:
            best = t
    return best

seen, reps = set(), []
for bits in range(1 << len(pairs)):
    A = np.zeros((n, n), dtype=np.int8)
    for k, (i, j) in enumerate(pairs):
        if bits >> k & 1: A[i, j] = 1
        else:             A[j, i] = 1
    if not (A.sum(axis=1) == 3).all():
        continue
    c = canon(A)
    if c in seen:
        continue
    seen.add(c)
    reps.append(A.copy())

print("regular tournaments on 7 vertices, up to isomorphism: %d" % len(reps))
print()
for idx, A in enumerate(reps, 1):
    M = sp.Matrix(A.tolist())
    lam = sp.Symbol('lam')
    cpA = sp.expand(M.charpoly(lam).as_expr())
    # substitute lambda = (x - 1)/2 and clear the 2^n
    cpS = sp.expand(sp.simplify(cpA.subs(lam, (x - 1) / 2) * 2**n))
    print("  [%d] char_A(lam) = %s" % (idx, sp.factor(cpA)))
    print("      char_S(x)   = %s" % sp.factor(cpS))
    fl = sp.factor_list(cpS)
    print("      factors     = %s" % [sp.expand(f) for f, _ in fl[1]])
    # is it x * (even) ?
    q = sp.simplify(sp.cancel(cpS / x))
    even = sp.simplify(sp.expand(q.subs(x, -x) - q)) == 0
    print("      char_S = x * EVEN polynomial : %s" % even)
    print()
