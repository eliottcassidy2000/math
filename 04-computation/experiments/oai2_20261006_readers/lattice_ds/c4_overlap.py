# FINITE-EXACT: pair-correlation of the LRC bad sets B_v = {t : ||v t|| < delta} (unrestricted numerators,
# exactly the sets E_v of #022 with gamma = 0, psi(v) = delta).  |B_v| = 2 delta.
# |B_q cap B_r| = d * sum_k g(k/L), d = gcd, L = lcm, g = overlap of [-delta/q,delta/q] and [x-delta/r,x+delta/r].
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def pair(q, r, dl):
    d = gcd(q, r); L = q*r//d; A = dl/q; B = dl/r
    tot = F(0); k = 0
    # sum over k in Z with |k/L| < A + B
    kmax = int((A + B)*L) + 1
    for k in range(-kmax, kmax+1):
        x = F(k, L)
        ov = min(A, x + B) - max(-A, x - B)
        if ov > 0: tot += ov
    return d*tot
def union_measure(S, dl, grid=None):
    # exact: complement of union = {t: ||v t|| >= dl for all v}; compute via breakpoints
    pts = {F(0), F(1)}
    for v in S:
        for a in range(v+1):
            for s in (-1, 1):
                x = (F(a) + s*dl)/v
                if 0 <= x <= 1: pts.add(x)
    pts = sorted(pts); good = F(0)
    for lo, hi in zip(pts, pts[1:]):
        mid = (lo + hi)/2
        if all(min((v*mid) % 1, 1 - (v*mid) % 1) >= dl for v in S): good += hi - lo
    return 1 - good
def report(name, S, dl):
    n = len(S); sumB = 2*dl*n
    pairs = sum(pair(q, r, dl) for q, r in combinations(S, 2))
    indep_pairs = sum((2*dl)**2 for _ in combinations(S, 2))
    U = union_measure(S, dl)
    indepU = 1 - (1 - 2*dl)**n
    print(f"{name:28s} n={n:2d} delta={str(dl):7s} sum|B|={float(sumB):.4f} sum pair-overlaps={float(pairs):.4f} "
          f"(indep {float(indep_pairs):.4f}, ratio {float(pairs/indep_pairs):.3f})  |union|={float(U):.4f} (indep model {float(indepU):.4f})")
dl = F(14, 183)
report("deep well {1..12,182}", list(range(1, 13)) + [182], dl)
report("AP {1..13}", list(range(1, 14)), F(1, 14))
report("AP {1..13} at 14/183", list(range(1, 14)), dl)
report("{1..11,13,84}", list(range(1, 12)) + [13, 84], dl)
random.seed(1)
for trial in range(3):
    S = sorted(random.sample(range(1, 200), 13))
    report(f"random 13-set #{trial}", S, dl)
