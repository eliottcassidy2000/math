#!/usr/bin/env python3
"""Rotation lemma: for a cyclic arrangement c of [n] and the identity linear order, the n rotations (cut points) of the
cycle give paths; the average over rotations and over pairs i<j of [i before j] is (1/(n*C(n,2))) * sum_{i<j} d(j -> i),
d = forward cyclic distance.  Claim: the maximum over all cyclic arrangements is attained by the identity cycle and equals
(2n-1)/(3n) < 2/3.  Brute force over all (n-1)! cyclic arrangements, n <= 9."""
from itertools import permutations
from fractions import Fraction as F
for n in range(3, 10):
    best = None; arg = None
    for perm in permutations(range(1, n)):
        cyc = (0,) + perm
        pos = {v: k for k, v in enumerate(cyc)}
        s = sum((pos[i] - pos[j]) % n for i in range(n) for j in range(i + 1, n))   # d(j -> i) forward from j to i
        if best is None or s > best: best, arg = s, cyc
    avg = F(best, n * (n * (n - 1) // 2))
    print(f"n={n}: max average {avg} = {float(avg):.4f}; (2n-1)/(3n) = {F(2*n-1, 3*n)}; argmax cycle {arg}")
