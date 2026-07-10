# kind-pasteur-2026-07-10-S127
# Anatomy of the TIGHT LOCUS (mu = 0) -- what actually excludes it from the residual class.
#
# mu(S) = meas {t in [0,1) : ||v*t|| >= 1/14 for all v in S}, computed EXACTLY (Fraction breakpoint sweep).
# Breakpoints for speed v: t with frac(v t) in {1/14, 13/14}  =>  t = (14k+1)/(14v), (14k+13)/(14v).
#
# FINDING (corrects the loose claim "covering => mu > 0"):
#   the dilate 2*{1..13} is COVERING and has mu = 0.
#   What excludes the tight locus is the SCALE GAP (GapFamily: ratio > 13) / primitivity -- not covering.
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

def safe_measure(S):
    """exact Lebesgue measure of the safe set, and its number of components"""
    bps = {F(0), F(1)}
    for v in S:
        for k in range(0, v):
            bps.add(F(14 * k + 1, 14 * v)); bps.add(F(14 * k + 13, 14 * v))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0); ncomp = 0; inprev = False
    for a, b in zip(bps, bps[1:]):
        mid = (a + b) / 2
        ok = all(min((v * mid) % 1, 1 - (v * mid) % 1) >= F(1, 14) for v in S)
        if ok: tot += (b - a); ncomp += (0 if inprev else 1)
        inprev = ok
    return tot, ncomp

def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S) == 1
def ratio_gt_13(S): return max(S) > 13 * min(S)      # GapFamily

print("=== THE TIGHT LOCUS (mu = 0): the dilated APs c*{1..13} ===")
for c in [1, 2, 3, 5]:
    S = [c * i for i in range(1, 14)]
    m, k = safe_measure(S)
    print(f"  {c}*{{1..13}}: mu={m}  components={k}  covering={covering(S)}  "
          f"primitive={primitive(S)}  ratio>13={ratio_gt_13(S)}  (max/min = {max(S)//min(S)})")

print()
print("=== THE CORRECTION: covering does NOT exclude mu = 0 ===")
S = [2 * i for i in range(1, 14)]
print(f"  2*{{1..13}} = {S}")
print(f"    covering   : {covering(S)}   <- YES")
print(f"    mu         : {safe_measure(S)[0]}   <- ZERO")
print(f"    primitive  : {primitive(S)}   <- NO  (gcd = 2)")
print(f"    ratio > 13 : {ratio_gt_13(S)}   <- NO  (max/min = 13 exactly)")
print("    => it is dispatched by spread13 (NOT GapFamily), never reaching the residual.")
print("  So 'covering => mu > 0' is FALSE. What excludes the tight locus is the SCALE GAP / primitivity.")

print()
print("=== every tight family has ratio EXACTLY 13 (so GapFamily excludes them all) ===")
for c in [1, 2, 3, 7]:
    S = [c * i for i in range(1, 14)]
    print(f"  {c}*{{1..13}}: max/min = {max(S)}/{min(S)} = {max(S)//min(S)}  (GapFamily needs > 13)")

print()
print("=== how small does mu get on covering families? (per-family floor, not uniform) ===")
import random; random.seed(2)
fams = [list(S) for S in combinations(range(1, 23), 13) if covering(list(S))]
vals = sorted((safe_measure(S)[0], S) for S in random.sample(fams, 30))
print(f"  sampled 30 of {len(fams)} covering 13-subsets of [1,22]; zero-mu count: "
      f"{sum(1 for m, _ in vals if m == 0)}")
for m, S in vals[:5]:
    print(f"    mu = {float(m):.6f} = {m}   {S}")
print("  (klein THM-685 Cor 1: any mu > 0 gives live rulers at ALL q > Sum(v)/mu -- the bank below is finite.")
print("   So NO uniform floor is needed; a per-family mu > 0 suffices.)")
