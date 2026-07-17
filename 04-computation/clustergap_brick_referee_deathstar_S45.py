#!/usr/bin/env python3
"""Referee for LRCClusterGapBrick.lean (death-star-S45): the corrected L1.

1. sorted_gap_pigeonhole: random rational teeth over random [a,b]; whenever
   0 < L - B (B = total clipped mass), a closed free segment of width
   >= (L-B)/(N+1) avoiding every open tooth must exist.
2. The counterexample: [0,1], tooth (-1,2), B=1 -> NO avoiding pair.
"""
from fractions import Fraction as F
import random
random.seed(45)

def free_segments(a, b, teeth):
    """UNMERGED breakpoint segments of [a,b] avoiding every open tooth in the
    ENDPOINT sense (hi <= t1 or t2 <= lo) -- the exact shape of the Lean
    conclusion.  No merging: merging across an EMPTY tooth (t1 == t2) builds a
    segment that straddles its endpoint and falsely fails the endpoint form;
    each tooth (empty or not) adds at most one cut, so the (N+1)-pigeonhole
    still governs the unmerged widths."""
    pts = sorted({a, b} | {max(min(x, b), a) for t in teeth for x in t})
    return [(lo, hi) for lo, hi in zip(pts, pts[1:])
            if not any(t1 < (lo + hi) / 2 < t2 for t1, t2 in teeth)]

trials, checked, fails = 6000, 0, 0
for _ in range(trials):
    a = F(random.randint(-20, 20), random.randint(1, 9))
    b = a + F(random.randint(1, 40), random.randint(1, 9))
    N = random.randint(0, 7)
    teeth = []
    for _ in range(N):
        t1 = a + F(random.randint(-15, 60), 12) * (b - a) / 4
        t2 = t1 + F(random.randint(0, 30), 17)
        teeth.append((t1, t2))
    B = sum(max(F(0), min(t2, b) - max(t1, a)) for t1, t2 in teeth)
    L = b - a
    if L - B <= 0:
        continue
    checked += 1
    need = (L - B) / (N + 1)
    segs = free_segments(a, b, teeth)
    best = max((hi - lo for lo, hi in segs), default=F(0))
    # verify the found best segment truly avoids every tooth (endpoint check)
    ok_avoid = True
    for lo, hi in segs:
        for t1, t2 in teeth:
            if not (hi <= t1 or t2 <= lo):
                ok_avoid = False
    if best < need or not ok_avoid:
        fails += 1
        if fails <= 3:
            print("FAIL", a, b, teeth, "best", best, "need", need)
print(f"sorted_gap_pigeonhole: {checked} positive-slack trials of {trials}, fails={fails}",
      "PASS" if fails == 0 else "FAIL")

# counterexample: [0,1], tooth (-1,2): every point of [0,1] is inside the tooth
cex_covered = all(not (F(x,100) <= -1 or 2 <= F(x,100)) for x in range(0, 101))
print("positivity-necessity counterexample: tooth (-1,2) covers all of [0,1]:",
      "PASS" if cex_covered else "FAIL")
