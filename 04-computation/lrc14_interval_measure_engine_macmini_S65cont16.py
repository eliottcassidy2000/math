#!/usr/bin/env python3
"""
lrc14_interval_measure_engine_macmini_S65cont16.py -- THE FINITE INTERVAL-MEASURE ENGINE
for the hsmall/hB legs: slowmu(safeSet P) EXACT for every P subset of {1..13}.

safeSet P = {x : forall p in P, fract(p x) in [1/14, 13/14]}  (LRCDenseCovers, verbatim);
per-runner: safe(p) cap [0,1) = Union_j [(14j+1)/(14p), (14j+13)/(14p)] -- rational intervals;
intersections stay rational-interval unions; measure = sum of lengths. Pure Fractions.

VALIDATION: reproduces THM-530's m_P = 14249/252252 at the exact argmin
{1,2,3,5,7,8,9,11,12,13} (min over |P| <= 10), from scratch.

FIRST OUTPUT = A SOUNDNESS FINDING (flagged, not overridden): the Lean skeleton's hsmall --
and hence the top-level hfloor (witnessMP <= witnessG2 for ALL v) -- is UNSATISFIABLE as
stated: v = (1,...,13) has clusterSize(shapeOf v) = 0 and witnessG2 = meas(safeSet{1..13}) = 0
< witnessMP. Exact failure boundary: |P| = 11: min 313/9702; |P| = 12: 7/858; |P| = 13: 0 --
i.e. the m_P floor fails EXACTLY for clusterSize <= 2, matching THM-530's admissibility
(k >= 3) which the skeleton statement dropped. REPAIR (proposed): hfloor needs only
POSITIVITY per shape to feed hpartA: k >= 3 -> m_P; k = 2 -> 313/9702; k = 1 -> 7/858;
k = 0 -> the single family v = +-{1..13} = the AP = non-covering -> the q = 14 sieve
(LonelyRunner.sieve_one_div). All constants exact, native_decide-shaped.
"""
from fractions import Fraction as F
from itertools import combinations
def bands(p): return [(F(14*j+1,14*p), F(14*j+13,14*p)) for j in range(p)]
def inter(A,B):
    out=[]
    for a1,b1 in A:
        for a2,b2 in B:
            lo,hi=max(a1,a2),min(b1,b2)
            if lo<hi: out.append((lo,hi))
    return sorted(out)
def measure(P):
    cur=[(F(0),F(1))]
    for p in sorted(P):
        cur=inter(cur,bands(p))
        if not cur: return F(0)
    return sum(b-a for a,b in cur)
if __name__ == "__main__":
    mP=F(14249,252252)
    best=None
    for size in range(1,11):
        for P in combinations(range(1,14),size):
            m=measure(P)
            if best is None or m<best[0]: best=(m,P)
    print(f"min |P|<=10: {best[0]} at {best[1]} (m_P reproduced: {best[0]==mP})")
    for size in range(1,14):
        mn=min((measure(P),P) for P in combinations(range(1,14),size))
        print(f"|P|={size:2d}: min = {str(mn[0]):>14} = {float(mn[0]):.6f} at {mn[1]}")
