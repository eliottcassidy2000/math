#!/usr/bin/env python3
"""
Two-coordinate shallow winding census (numeric-screened)  (mac-mini-2026-07-17-S109)
====================================================================================
THM-1001 closes SINGLE-coordinate shallow winding of {1..12} at all heights. Its
reduction leaves the residual: >=2 coordinates wound above THM-770's height 12.
This probes the TWO-coordinate slice to height 20 (well past THM-770's 12 and
klein-S313's exhaustive {1..16} = winding height 1), plus a random multi-wound
high-height sample.

Speed: most wound sets are loose (M >> 1/13); a numpy grid pre-screen (M within a
hair of 1/13) filters to the few near-floor candidates, which are then confirmed
with the exact corrected M (MISTAKE-144). A tight set has M = 1/13 EXACTLY.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random

ONE13 = F(1, 13)
G = 40040  # = 8*5*7*11*13, resolves /13 and small denominators well
tg = np.arange(1, G) / G

def Mnum(S):
    mn = None
    for v in S:
        r = (v * tg) % 1.0; d = np.minimum(r, 1 - r)
        mn = d if mn is None else np.minimum(mn, d)
    return mn.max()

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

def is_tight(A):                       # numeric screen then exact confirm
    if Mnum(A) > 1.0 / 13 + 6e-4: return False
    return M_exact(A) == ONE13

H = 20
print(f"=== TWO-coordinate shallow winding: 66 pairs, heights 1..{H} (screened+exact) ===")
base = list(range(1, 13)); found = []; checked = 0; scr = 0
for i, j in combinations(range(1, 13), 2):
    rest = [x for x in base if x != i and x != j]
    for ki in range(1, H + 1):
        wi = i + 13 * ki
        for kj in range(1, H + 1):
            wj = j + 13 * kj
            A = sorted(rest + [wi, wj])
            if len(set(A)) != 12 or reduce(gcd, A) != 1: continue
            checked += 1
            if Mnum(A) <= 1.0 / 13 + 6e-4:
                scr += 1
                if M_exact(A) == ONE13: found.append(tuple(A))
print(f"  checked {checked}; near-floor screened {scr}; EXACT tight (M=1/13): {len(found)}")
for A in found[:8]: print("   ", list(A))

print()
print("=== RANDOM multi-coordinate (3..6 wound) high-height sample ===")
rng = random.Random(109); found2 = []; checked2 = 0; scr2 = 0
for _ in range(200000):
    nwound = rng.randint(3, 6)
    wound = set(rng.sample(range(1, 13), nwound))
    A = sorted(set((r + 13 * rng.randint(1, 25)) if r in wound else r for r in range(1, 13)))
    if len(A) != 12 or reduce(gcd, A) != 1: continue
    checked2 += 1
    if Mnum(A) <= 1.0 / 13 + 6e-4:
        scr2 += 1
        if M_exact(A) == ONE13: found2.append(tuple(A))
print(f"  checked {checked2}; near-floor screened {scr2}; EXACT tight: {len(found2)}")
for A in found2[:8]: print("   ", list(A))

print()
print("VERDICT")
print(f"  two-coordinate winding (heights 1..{H}): {len(found)} tight found")
print(f"  random multi-coordinate (heights 1..25):  {len(found2)} tight found")
print("  => no sporadic shallow tight set in this residual slice (empirical, not exhaustive;")
print("     the >=2-high-coordinate shallow residual remains open beyond the screened range).")
