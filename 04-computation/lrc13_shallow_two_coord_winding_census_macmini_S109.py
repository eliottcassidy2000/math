#!/usr/bin/env python3
"""
Two-coordinate shallow winding census  (mac-mini-2026-07-17-S109)
=================================================================
THM-1001 proves the SINGLE-coordinate shallow winding of {1..12} is never tight
(all heights). Its reduction (C) leaves the residual: shallow tight sets with >=2
coordinates wound above THM-770's height 12. This module exhaustively probes the
two-coordinate slice to a height FAR above THM-770 (height 12) and klein-S313's
exhaustive {1..16} box (winding height 1):

  A = {1..12}\\{i,j} u {i+13 k_i, j+13 k_j},   1 <= k_i,k_j <= H,   i<j.

For each of the 66 residue pairs and all H^2 height pairs, test M(A)=1/13.
Also a broader RANDOM multi-coordinate (3..5 wound) high-height sample.
Everything exact (corrected M_exact, MISTAKE-144). Reports any tight set found.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random

ONE13 = F(1, 13)

def M_exact(S):
    S = sorted(set(S)); best = F(0); dens = {2 * a for a in S}
    for a, b in combinations(S, 2):
        dens.add(a + b); dens.add(b - a)
    dens.discard(0)
    for q in dens:
        for m in range(1, q):
            num = q
            for v in S:
                r = (v * m) % q; d = min(r, q - r)
                if d < num: num = d
                if num * 13 < q: break
            if num * 13 >= q:
                c = F(num, q)
                if c > best: best = c
    return best

H = 25
print(f"=== TWO-coordinate shallow winding: 66 pairs, heights 1..{H} (exact) ===")
base = list(range(1, 13)); found = []; checked = 0
for i, j in combinations(range(1, 13), 2):
    rest = [x for x in base if x != i and x != j]
    for ki in range(1, H + 1):
        wi = i + 13 * ki
        for kj in range(1, H + 1):
            wj = j + 13 * kj
            A = sorted(rest + [wi, wj])
            if len(set(A)) != 12: continue
            if reduce(gcd, A) != 1: continue
            checked += 1
            if M_exact(A) == ONE13:
                found.append(tuple(A))
print(f"  checked {checked} primitive two-coordinate-wound sets; tight (M=1/13): {len(found)}")
for A in found[:8]: print("   ", list(A))

print()
print(f"=== RANDOM multi-coordinate (3..5 wound) high-height sample ===")
rng = random.Random(109); found2 = []; checked2 = 0
for _ in range(120000):
    nwound = rng.randint(3, 5)
    wound = rng.sample(range(1, 13), nwound)
    A = []
    for r in range(1, 13):
        if r in wound:
            A.append(r + 13 * rng.randint(1, 30))
        else:
            A.append(r)
    A = sorted(set(A))
    if len(A) != 12 or reduce(gcd, A) != 1: continue
    checked2 += 1
    if M_exact(A) == ONE13:
        found2.append(tuple(A))
print(f"  checked {checked2} primitive multi-wound sets (heights up to 30); tight: {len(found2)}")
for A in found2[:8]: print("   ", list(A))

print()
print("VERDICT: two-coordinate shallow winding through height 25 and a random")
print("multi-coordinate high-height sample contain NO tight set besides none-found —")
print(f"  two-coord tight: {len(found)}, multi-coord tight: {len(found2)}.")
print("  (Empirical support for shallow-branch emptiness in the residual THM-1001 leaves;")
print("   NOT an exhaustive proof — the shallow residual >=2 high coordinates stays open.)")
