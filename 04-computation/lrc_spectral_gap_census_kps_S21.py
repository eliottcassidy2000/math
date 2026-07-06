#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S21: BROAD SPECTRAL-GAP CENSUS -- independent verification
that no 12-speed family has M(v) in the OPEN gap (1/13, 2/25).

The whole LRC(14) crux (after the fleet's reductions: J-K / mac-mini subsumption +
opus residue bridge) collapses to ONE finite object: the 1-D census -- no primitive
12-family v has M(v) = max_t min_i ||v_i t|| strictly inside (1/13, 2/25).
  1/13 = 0.076923...   2/25 = 0.080000

mac-mini ran a STRUCTURED template census (tight families c*{1..12} etc.).  This is
the COMPLEMENTARY brute force: ALL 12-subsets of speeds up to a height bound, EXACT
rational M, checking none lands in the open gap.  Independent lane, same object.

EXACT M: M(v) = max over rationals a/q of min_i ||v_i a/q||
       = max_{q<=Q} max_{a} [ min_i min(v_i a mod q, q - v_i a mod q) ] / q.
For lonely-runner-type max-min the optimum sits at a bounded denominator; we take
Q = 4*sum|v_i| (validated on known families) and report exact fractions.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce
import numpy as np

LO = Fraction(1, 13)   # 0.076923...
HI = Fraction(2, 25)   # 0.080000

def tuple_gcd(v):
    return reduce(gcd, (abs(x) for x in v))

def M_exact(v, Qcap=None):
    """Exact M(v) = max_{q,a} min_i ||v_i a/q||.  Returns Fraction (numpy-vectorized)."""
    S = int(sum(abs(x) for x in v))
    Q = Qcap if Qcap else 4 * S
    va = np.array(v, dtype=np.int64)
    best_num, best_den = 0, 1   # best fraction = best_num/best_den
    for q in range(2, Q + 1):
        a = np.arange(1, q, dtype=np.int64)              # (q-1,)
        r = np.outer(va, a) % q                          # (l, q-1)
        d = np.minimum(r, q - r)                         # (l, q-1)
        bestq = int(d.min(axis=0).max())                 # max_a min_i
        # compare bestq/q > best_num/best_den
        if bestq * best_den > best_num * q:
            best_num, best_den = bestq, q
    return Fraction(best_num, best_den)

# ---- validation on KNOWN critical families ----
print("=== validation on known families (12 speeds) ===", flush=True)
known = {
    "AP {1..12}": list(range(1, 13)),
    "{1..11,24}": list(range(1, 12)) + [24],
    "{1..11,13}": list(range(1, 12)) + [13],
    "geometric-ish {1,2,3,5,7,8,9,11,13,17,19,23}": [1,2,3,5,7,8,9,11,13,17,19,23],
}
for name, v in known.items():
    M = M_exact(v)
    ingap = LO < M < HI
    print(f"  {name}: M = {M} = {float(M):.6f}  {'*** IN GAP ***' if ingap else '(ok: '+('<=1/13' if M<=LO else '>=2/25')+')'}", flush=True)

# ---- broad brute force: all 12-subsets of {1..H}, primitive, check gap ----
print(flush=True)
print("=== broad census: 12-subsets of {1..H}, primitive, M in open gap? ===", flush=True)
import sys
gap_hits = []
near = []  # M within [0.074, 0.082] -- near the gap, worth reporting
for H in [13, 14, 15, 16, 17, 18]:
    count = 0
    minM_above = Fraction(1)   # smallest M that is >= 2/25
    maxM_below = Fraction(0)   # largest M that is <= 1/13
    for combo in combinations(range(1, H + 1), 12):
        v = list(combo)
        if tuple_gcd(v) != 1:
            continue
        count += 1
        M = M_exact(v)
        if LO < M < HI:
            gap_hits.append((v, M))
        elif M >= HI and M < minM_above:
            minM_above = M
        elif M <= LO and M > maxM_below:
            maxM_below = M
        if Fraction(74,1000) <= M <= Fraction(82,1000) and not (LO < M < HI):
            near.append((v, M))
    print(f"  H={H}: {count} primitive 12-subsets; gap hits so far = {len(gap_hits)}; "
          f"tightest below 1/13 = {float(maxM_below):.6f}, tightest above 2/25 = {float(minM_above):.6f}", flush=True)

print(flush=True)
if gap_hits:
    print(f"*** {len(gap_hits)} GAP MEMBERS FOUND (would be a counterexample!) ***", flush=True)
    for v, M in gap_hits[:10]:
        print(f"    v={v}  M={M}={float(M):.6f}", flush=True)
else:
    print("NO gap members in the broad census -- the open gap (1/13, 2/25) is EMPTY over all", flush=True)
    print("primitive 12-subsets of {1..18}.  Independent confirmation of the census object.", flush=True)
print(flush=True)
print(f"near-gap families (M in [0.074, 0.082], on either edge): {len(near)}", flush=True)
for v, M in sorted(near, key=lambda x: abs(float(x[1]) - 0.078))[:12]:
    edge = "1/13 edge" if M <= LO else "2/25 edge"
    print(f"    v={v}  M={float(M):.6f}  ({edge})", flush=True)
