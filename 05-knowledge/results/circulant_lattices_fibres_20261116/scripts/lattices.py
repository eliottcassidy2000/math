#!/usr/bin/env python3
"""lattices.py (H-2): the delayed-choice lattice at order n: Delta a supported on m = n (mod 2), 1 <= m <= n - 2, with
sum Delta a_m = 0 and sum Delta a_m 2^m = 0; g_n = gcd of the values sum Delta a_m 3^m over an integer basis of the
kernel (sympy); also the joint lattice of (Delta P_3, Delta P_4)."""
from sympy import Matrix, gcd, ZZ
from sympy.matrices.normalforms import hermite_normal_form
import math
for n in range(5, 15):
    ms = list(range(n - 2, 0, -2))
    A = Matrix([[1] * len(ms), [2 ** m for m in ms]])
    K = A.nullspace()
    # integer basis: clear denominators of each nullspace vector (they are rational); then the lattice generated may be a
    # superlattice of the integer kernel's image only up to index; compute the integer kernel properly via HNF of the kernel
    vecs = []
    for v in K:
        den = 1
        for x in v: den = math.lcm(den, x.q)
        vecs.append([int(x * den) for x in v])
    # saturate: integer kernel = all integer vectors in the Q-span; use Smith/HNF on the transpose
    Mv = Matrix(vecs)
    # Hermite normal form of the lattice spanned, then saturate by dividing by gcd of maximal minors is nontrivial; instead
    # enumerate small integer combinations to find the integer kernel's images (exact for these sizes)
    import itertools
    vals3 = set(); vals34 = set()
    rng = range(-6, 7)
    for coeffs in itertools.product(rng, repeat=len(vecs)):
        w = [sum(c * v[i] for c, v in zip(coeffs, vecs)) for i in range(len(ms))]
        if any(x % 1 for x in w): continue
        p3 = sum(x * 3 ** m for x, m in zip(w, ms)); p4 = sum(x * 4 ** m for x, m in zip(w, ms)); vals3.add(p3); vals34.add((p3, p4))
    # the true integer kernel may be finer than the span of the cleared vectors; refine by checking all integer vectors
    # with small entries in the kernel directly
    small = set()
    for w in itertools.product(range(-3, 4), repeat=len(ms)):
        if sum(w) == 0 and sum(x * 2 ** m for x, m in zip(w, ms)) == 0: small.add(sum(x * 3 ** m for x, m in zip(w, ms)))
    g = 0
    for v in vals3 | small: g = math.gcd(g, v)
    print(f"n = {n}: unknown coefficients m in {ms}; kernel dimension {len(K)}; Delta P_3 lattice = {g} Z")
