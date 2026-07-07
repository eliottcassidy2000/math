"""
mac-mini-2026-07-07-S52 (HYP-5107, THM-649) -- THE CYLINDER PARITY LAW (the pair-space
dimension count) + verification for the canon trilogy fold.

DERIVATION (to be verified before assertion):
  Crossing count between two winding cross-edges: N = #(Z cap open(x, x+delta)),
  delta = const_pair + (w1 - w2).  For generic twist (x, x+delta never integers),
  shifting any single w by 1 moves one endpoint of the interval by 1 => the count
  changes by EXACTLY 1 => mod 2, every pair contributes c0_pair + w1 + w2.
  Summing over pairs:  Q_cross(w) == Q_cross(0) + (E-1) * sum_e w_e  (mod 2),
  E = #cross edges = m*n'.  THE PAIR-SPACE (quadratic part mod 2) HAS DIMENSION 0.
  Corollaries:
    - E odd  => (E-1) even => CROSSING PARITY INVARIANT over the whole winding cube
      (all cylindrical drawings with these spines share Q mod 2) -- for K_{m,n} with
      m, n odd this is the cylindrical Kleitman parity phenomenon;
    - E even => Q mod 2 = const + total winding parity (a single linear functional).

VERIFY: exhaustive over winding cubes at several (m, n', twist): the mod-2 residuals
of Q against the affine model must vanish; report the invariant/functional per case.
"""
import numpy as np
from itertools import combinations
import math, random as rnd
rnd.seed(52)

def cross_count(a1, d1, a2, d2):
    x = a1 - a2; delta = d1 - d2
    lo, hi = min(x, x+delta), max(x, x+delta)
    if hi <= lo: return 0
    return math.floor(hi) - math.floor(lo) - (1 if hi == math.floor(hi) else 0)

def Qcross(m, np_, wflat, tw):
    edges = [(i, j) for i in range(m) for j in range(np_)]
    tot = 0
    for (e1, (i1, j1)), (e2, (i2, j2)) in combinations(list(enumerate(edges)), 2):
        a1, b1 = i1/m, (j1 + tw)/np_
        a2, b2 = i2/m, (j2 + tw)/np_
        d1 = b1 + (-wflat[e1]) - a1     # w in {0,1} -> winding {0,-1}
        d2 = b2 + (-wflat[e2]) - a2
        tot += cross_count(a1, d1, a2, d2)
    return tot

print("=== the cylinder parity law: exhaustive mod-2 verification ===")
for (m, np_) in [(2,3), (3,3), (2,2), (3,4), (2,4)]:
    E = m*np_
    tw = 0.137  # generic twist
    base = Qcross(m, np_, [0]*E, tw)
    ok_affine = True; inv = True
    for mask in range(1 << E):
        w = [(mask >> b) & 1 for b in range(E)]
        q = Qcross(m, np_, w, tw)
        pred = (base + (E-1)*sum(w)) % 2
        if q % 2 != pred: ok_affine = False
        if q % 2 != base % 2: inv = False
    print(f"  (m,n')=({m},{np_}) E={E} ({'odd' if E%2 else 'even'}): affine law {'HOLDS' if ok_affine else 'FAILS'}"
          f" on all 2^{E}; parity {'INVARIANT' if inv else 'varies (= total-winding functional)' }"
          f"  [prediction: {'invariant' if (E-1)%2==0 else 'functional'}]")

# random-twist robustness
print("\n  twist robustness (m,n')=(3,3), 5 random twists:")
for _ in range(5):
    tw = rnd.random()
    base = Qcross(3, 3, [0]*9, tw)
    bad = 0
    for mask in range(1 << 9):
        w = [(mask >> b) & 1 for b in range(9)]
        if Qcross(3, 3, w, tw) % 2 != base % 2: bad += 1
    print(f"    tw={tw:.3f}: parity-invariance violations {bad}/512")
