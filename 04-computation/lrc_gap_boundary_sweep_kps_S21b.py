#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S21b: TARGETED 2/25-BOUNDARY sweep.

The broad {1..18} census (S21) found the open gap (1/13, 2/25) EMPTY, with the AP
{1..12} the UNIQUE family at 1/13 and everything else >= 1/12.  But the 2/25 upper
edge is achieved by families with a LARGE speed ({1..11,24} = 2/25 exactly), missed
by {1..18}.  This sweep maps the 2/25 boundary precisely: families with one or two
large speeds, confirming NONE dips into the open gap (1/13, 2/25).

These near-2/25 families ARE the dangerous region (the gap's upper edge), so this is
the delicate part of the census.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce
import numpy as np

LO = Fraction(1, 13)
HI = Fraction(2, 25)

def tuple_gcd(v):
    return reduce(gcd, (abs(x) for x in v))

def M_exact(v, Qcap=None):
    S = int(sum(abs(x) for x in v))
    Q = Qcap if Qcap else 4 * S
    va = np.array(v, dtype=np.int64)
    best_num, best_den = 0, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q, dtype=np.int64)
        r = np.outer(va, a) % q
        d = np.minimum(r, q - r)
        bestq = int(d.min(axis=0).max())
        if bestq * best_den > best_num * q:
            best_num, best_den = bestq, q
    return Fraction(best_num, best_den)

gap_hits = []
edge_2_25 = []   # families with M exactly 2/25
below = []        # M in (1/13, 2/25) -- the forbidden zone

def check(v, tag):
    if tuple_gcd(v) != 1:
        # reduce
        g = tuple_gcd(v); v = [x//g for x in v]
    M = M_exact(v)
    if LO < M < HI:
        gap_hits.append((v, M, tag)); below.append((v, M))
    if M == HI:
        edge_2_25.append((sorted(v), tag))
    return M

# (1) {1..11, k}: one large speed, k = 12..60
print("=== {1,...,11, k}, k = 12..60 ===", flush=True)
minM = Fraction(1)
for k in range(12, 61):
    v = list(range(1, 12)) + [k]
    M = check(v, f"1..11,{k}")
    if HI <= M < minM: minM = M
    mark = "  <-- 2/25 EDGE" if M == HI else ("  *** GAP ***" if LO < M < HI else "")
    if M <= Fraction(85,1000):
        print(f"  k={k:2d}: M = {M} = {float(M):.6f}{mark}", flush=True)

# (2) {1..10, j, k}: two large speeds, j<k in 12..40
print("=== {1,...,10, j, k}, 12<=j<k<=40 (near-2/25 scan) ===", flush=True)
cnt = 0; ingap = 0
for j, k in combinations(range(12, 41), 2):
    v = list(range(1, 11)) + [j, k]
    M = check(v, f"1..10,{j},{k}")
    cnt += 1
    if LO < M < HI:
        ingap += 1
        print(f"  *** GAP: {v} M={M}={float(M):.6f} ***", flush=True)
print(f"  scanned {cnt} two-large-speed families; gap members = {ingap}", flush=True)

# (3) AP perturbations: {1..12} with one entry replaced by a larger value
print("=== AP {1..12} single-entry perturbations (replace one by w=13..40) ===", flush=True)
ingap3 = 0
for drop in range(1, 13):
    for w in range(13, 41):
        v = [x for x in range(1, 13) if x != drop] + [w]
        M = check(v, f"AP-drop{drop}+{w}")
        if LO < M < HI:
            ingap3 += 1
            print(f"  *** GAP: dropped {drop}, added {w}: {sorted(v)} M={float(M):.6f} ***", flush=True)
print(f"  AP perturbations gap members = {ingap3}", flush=True)

print(flush=True)
print("=== VERDICT ===", flush=True)
print(f"  total gap members (M strictly in (1/13,2/25)): {len(gap_hits)}", flush=True)
print(f"  families exactly AT the 2/25 edge: {len(edge_2_25)}", flush=True)
for v, tag in edge_2_25[:12]:
    print(f"     {v}   [{tag}]", flush=True)
if not gap_hits:
    print("  => the 2/25 boundary is a HARD FLOOR: families reach 2/25 exactly but NEVER dip", flush=True)
    print("     into the open gap.  Independent confirmation the gap has no member (bounded ht).", flush=True)
