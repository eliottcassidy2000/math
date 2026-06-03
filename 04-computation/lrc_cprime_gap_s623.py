#!/usr/bin/env python3
"""
S623 — the QUANTITATIVE C'(n): min M(S) over primitive multiple-of-n configs.
C'(n) says M>1/n; the *quantitative* form asks HOW MUCH more — the second-loneliness gap
(Kravitz normalization).  We compute min M over an exhaustive box, watch it vs box size
(convergence), and compare to candidate clean values 1/(n-1), 2/(2n-1), s/(ns+1) ladder.
The argmin configs are the 'hardest' multiple-of-n sets — the real target for the
constructive route to LRC(n) (and the model for n=14).
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_n14_flowshells_s622 import gap_and_argmax

def min_M_multiple(n, M):
    lvl = Fr(1, n); best = None; arg = None
    for S in itertools.combinations(range(1, M+1), n-1):
        if not any(v % n == 0 for v in S): continue
        if reduce(gcd, S) != 1: continue
        m, _ = gap_and_argmax(S)
        if best is None or m < best:
            best, arg = m, S
    return best, arg

def kravitz_ladder(n, smax=4):
    return [Fr(s, n*s+1) for s in range(1, smax+1)]   # Kravitz-style s/(ns+1)

if __name__ == "__main__":
    print("min M(S) over primitive MULTIPLE-OF-n configs  (the quantitative C'(n) gap)")
    print("candidate clean values: 1/(n-1), 2/(2n-1); Kravitz ladder s/(ns+1)\n")
    for n in range(4, 9):
        boxes = [4*n, 6*n] if n <= 5 else ([3*n, 4*n] if n==6 else [3*n])
        row = []
        for M in boxes:
            g, arg = min_M_multiple(n, M)
            row.append((M, g, arg))
        last = row[-1]
        import sys as _s
        print(f" n={n}:  1/n={Fr(1,n)}={float(Fr(1,n)):.4f}  1/(n-1)={Fr(1,n-1)}  2/(2n-1)={Fr(2,2*n-1)}={float(Fr(2,2*n-1)):.4f}")
        for M, g, arg in row:
            print(f"    box<= {M:3d}:  min M = {g} = {float(g):.5f}   argmin {arg}")
        gv = last[1]
        # is it on the Kravitz ladder s/(ns+1)?
        s = gv/(1 - n*gv) if (1-n*gv) != 0 else None
        onlad = (s is not None and s.denominator == 1)
        print(f"    -> min M = {gv}; on Kravitz ladder s/(ns+1)? {onlad} (s={s});  ==2/(2n-1)? {gv==Fr(2,2*n-1)};  ==1/(n-1)? {gv==Fr(1,n-1)}")
