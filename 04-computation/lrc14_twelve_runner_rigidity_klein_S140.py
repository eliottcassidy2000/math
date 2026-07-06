#!/usr/bin/env python3
"""
klein-2026-07-05-S140 - THE 12-RUNNER RIGIDITY + SPECTRAL GAP (the loose branch, exact math first).

The open crux of LRC(14) is the LOOSE BRANCH of TightLooseDichotomy: every primitive compressed
covering family whose 12-runner BASE (the non-max runners) is NOT a dilated AP c*{1..12} has a real
witness t with min_i ||v_i t|| >= 2/25.  Equivalently (12-runner rigidity):

    M(B) := max_t min_{v in B} ||v t||  =  1/13   iff  B = c*{1..12}  (dilated AP);
                                          >= 2/25  otherwise.

i.e. the 12-runner covering-min spectrum has an EMPTY GAP (1/13, 2/25), with 1/13 attained ONLY by
dilated APs.  (2/25 = 0.08; 1/13 = 0.076923; the gap = 2/25 - 1/13 = 1/325.)

This script establishes the EXACT statement to prove: enumerate 12-runner families (dilated APs,
single/double perturbations of the AP, covering-bases, small subsets), compute M EXACTLY, and:
  (a) confirm the gap (1/13, 2/25) is empty;
  (b) confirm 1/13 is attained ONLY by dilated APs c*{1..12};
  (c) find the exact 2/25-attainers (the second-value extremizers) -> the rigidity's "first perturbation";
  (d) map the low spectrum (the values just above 1/13).
Exact (Fractions).  Save to 05-knowledge/results/.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

ONE13 = F(1, 13); TWO25 = F(2, 25)

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval(S, Qcap):
    best = F(0)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m
    return best

def tuple_gcd(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g

def qcap(S): return min(2 * max(S) + 4, 400)

print(f"1/13 = {float(ONE13):.6f}   2/25 = {float(TWO25):.6f}   gap = 2/25-1/13 = {TWO25-ONE13} = {float(TWO25-ONE13):.6f}")
print("=" * 84)

AP = list(range(1, 13))
# collect (M, family) for a broad set of 12-runner families
seen = {}
def add(S):
    S = tuple(sorted(S))
    if len(set(S)) != 12 or any(x <= 0 for x in S): return
    if S in seen: return
    seen[S] = Mval(list(S), qcap(list(S)))

# 1. dilated APs
for c in range(1, 9): add([c * x for x in AP])
# 2. single perturbations: drop one AP element, add x
for j in AP:
    for x in range(1, 40):
        if x not in AP or x == j:
            add([v for v in AP if v != j] + [x])
# 3. double perturbations: drop two, add two (small)
for drop in combinations(AP, 2):
    base = [v for v in AP if v not in drop]
    for x in range(13, 30):
        for y in range(x + 1, 31):
            add(base + [x, y])
# 4. all 12-subsets of [1,16] (near-complete small)
for S in combinations(range(1, 17), 12):
    add(list(S))
# 5. dilated near-APs: c*(AP with one perturbation)
for c in (2, 3):
    for j in AP:
        for x in range(1, 25):
            if x not in AP or x == j:
                add([c * v for v in ([w for w in AP if w != j] + [x])])

print(f"computed M for {len(seen)} distinct 12-runner families")
print("=" * 84)

# (a)+(b): who is in [1/13, 2/25)? and who attains exactly 1/13?
below = [(m, S) for S, m in seen.items() if ONE13 <= m < TWO25]
at13 = [(m, S) for S, m in seen.items() if m == ONE13]
print(f"families with M in [1/13, 2/25):  {len(below)}")
for m, S in sorted(below)[:20]:
    g = tuple_gcd(list(S)); red = tuple(x // g for x in S)
    isap = (red == tuple(AP))
    print(f"   M={str(m):>7} (~{float(m):.6f})  gcd={g}  {'DILATED-AP' if isap else 'NON-AP!!'}  {list(S)}")
print(f"\nfamilies attaining EXACTLY 1/13: {len(at13)}")
non_ap_at_13 = [(m,S) for m,S in at13 if tuple(x//tuple_gcd(list(S)) for x in S) != tuple(AP)]
print(f"   of which NON-dilated-AP: {len(non_ap_at_13)}  {'<-- RIGIDITY FAILS' if non_ap_at_13 else '(rigidity holds: 1/13 = dilated AP only)'}")
for m,S in non_ap_at_13[:8]: print(f"      {list(S)}")

# (c): the 2/25 attainers (second value)
at225 = [(m, S) for S, m in seen.items() if m == TWO25]
print(f"\nfamilies attaining EXACTLY 2/25 (the second value): {len(at225)}")
for m, S in sorted(at225)[:12]:
    g = tuple_gcd(list(S)); red = tuple(x // g for x in S)
    print(f"   {list(S)}   reduced={red}")

# (d): the low spectrum
vals = sorted(set(seen.values()))
print(f"\nLOW 12-runner spectrum (smallest M values found):")
for m in vals[:10]:
    ex = min([S for S, mm in seen.items() if mm == m])
    print(f"   M={str(m):>8} (~{float(m):.6f})   e.g. {list(ex)}")

print("\nREADING: if [1/13,2/25) contains ONLY dilated APs at exactly 1/13, the gap is confirmed and")
print("the loose-branch rigidity is: non-dilated-AP 12-base => M >= 2/25. The 2/25-attainers show the")
print("'first perturbation' structure (what the AP degrades into). That is the exact statement to prove.")
print("DONE")
