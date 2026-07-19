#!/usr/bin/env python3
"""
death-star-2026-07-19-S59 -- HYP-7885 part 3: the TWO-FAR NESTED ABSORPTION CLOSURE
of the N=12 first gap (CRUX (C)) Hamming-2 stratum.

Stratum: A = {1..12}\{i,j} u {x1,x2}, i<j, 13 <= x1 < x2 (the genuinely-two-new
case; x in {i,j} degenerates to the single-far stratum, closed by part 2).

NESTED ABSORPTION (theta = 2/25, the window top):
  L1 (two-far absorption): if x1 >= (1+2theta)/l(B) with B = {1..12}\{i,j}, then
     the longest base-safe interval contains a FULL inter-arc gap of the x1-comb
     (length (1-2theta)/x1), and since x2 > x1 > (2theta/(1-2theta)) x1 * (4/21<1)
     that gap contains an x2-safe point => M(A) >= theta.  So x1 < X1 :=
     ceil((1+2theta)/l(B)) is forced for any gap member.
  L2 (single-far absorption over the augmented base): for each x1 < X1, base
     B' = B u {x1} (11 speeds => M(B') >= 1/12 > theta by SETTLED LRC(<=13), so
     l(B') > 0): x2 >= X2 := ceil(2theta/l(B')) => M(A) >= theta.
  Remaining: x2 in (x1, X2) -- FINITE, checked exactly.
Every branch is decided; the stratum classification is COMPLETE.

Base floors: B has 10 speeds (M >= 1/11), B' has 11 (M >= 1/12); both > 2/25. OK.
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
from itertools import combinations
import sys, time

sys.path.insert(0, '04-computation')
from lrc_singlefar_absorption_atlas_deathstar_S59 import (
    M_exact, M_exact_wit, safe_intervals, cand_denoms)

def main():
    log = lambda s="": print(s, flush=True)
    N = 12
    theta = F(2, 25); lo = F(1, 13)
    log("== HYP-7885 part 3: two-far nested absorption closure, N=12 ==\n")
    members = []
    tot_checked = 0
    t0 = time.time()
    for i, j in combinations(range(1, N+1), 2):
        B = [v for v in range(1, N+1) if v not in (i, j)]
        iv, lB = safe_intervals(B, theta)
        X1 = ceil((1 + 2*theta) / lB)
        n_x1 = 0
        for x1 in range(13, X1):
            if x1 in B: continue
            n_x1 += 1
            B2 = B + [x1]
            iv2, lB2 = safe_intervals(B2, theta)
            if lB2 <= 0:
                log(f"  (i,j)=({i},{j}) x1={x1}: l(B')=0 UNEXPECTED"); continue
            X2 = ceil(2*theta / lB2)
            for x2 in range(x1+1, X2):
                if x2 in B2: continue
                m = M_exact(B2 + [x2], stop_above=theta)
                tot_checked += 1
                if lo < m < theta:
                    Mw, q, a = M_exact_wit(B2 + [x2])
                    members.append((i, j, x1, x2, Mw, q))
                    log(f"  !! MEMBER (i,j)=({i},{j}) x=({x1},{x2}) M={Mw} q={q}")
        log(f"  (i,j)=({i:>2},{j:>2}): l(B)={float(lB):.5f} X1={X1:>3} "
            f"x1-count={n_x1} cumulative x2-checks={tot_checked} "
            f"[{time.time()-t0:.0f}s]")
    log(f"\n== RESULT: two-far (Hamming-2) stratum of the N=12 gap ==")
    log(f"  exact checks: {tot_checked}; members: "
        f"{members if members else 'EMPTY -- stratum CLOSED'}")

if __name__ == "__main__":
    main()
