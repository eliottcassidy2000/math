# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont69: SUPPORT opus-S267's L2 large-sieve energy route -- pin the TRUE worst-case core*Sum eps^2.
#
# opus-S267: LRC(14)-covering <= (core * Sum_v eps_v^2) < (6/7)^2 = 0.735 via Cauchy-Schwarz (Sum|eps|<=sqrt(core*Sum eps^2)
# then coreCover<1). opus reports core*Sum eps^2 <= 0.328. But for |core|=1 (the residual, deep well etc.),
# core*Sum eps^2 = eps_1^2 = (coreCover - 1/7)^2, and cont.63 found coreCover up to ~0.92 => eps_1^2 ~ 0.60 >> 0.328.
# This pins the TRUE worst |core|=1 body and the ACTUAL margin to 0.735 that opus's tight large-sieve must deliver.
from math import gcd
from functools import reduce
from itertools import combinations
P = 2*3*5*7*11*13
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def is_cov(v,N=14): return all(any(x%d==0 for x in v) for d in range(2,N+1))
def prim(v): return reduce(gcd,v)==1
def core(v): return [x for x in v if gcd(x,P)==1]
def coreCover(B, level=1/14, grid=800000):
    """|D_1 cap G'|/|G'|: fraction of the body good-set G'={all b lonely} that lies in D_1={||t||<level}."""
    g=0; ins=0
    for k in range(grid):
        t=k/grid
        ok=True
        for b in B:
            r=(b*t)%1.0
            if min(r,1-r)<level: ok=False; break
        if ok:
            g+=1
            r0=t%1.0
            if min(r0,1-r0)<level: ins+=1
    return (ins/g) if g else 1.0

def main():
    tgt2 = (6/7)**2   # 0.7347 -- the L2 energy ceiling
    print(f"opus-S267 ceiling: core*Sum eps^2 < (6/7)^2 = {tgt2:.4f}; opus reports <= 0.328 (NOT the worst case).")
    print("For |core|=1: core*Sum eps^2 = eps_1^2 = (coreCover - 1/7)^2. Find the WORST |core|=1 covering body.\n")
    # candidate |core|=1 covering bodies B (12 non-core, {1}uB covering). Search structured tight ones.
    bodies = {}
    bodies["deep well {2..12,182}"]=list(range(2,13))+[182]
    bodies["{2..12,84}(84=6*14)"]=list(range(2,13))+[84]
    bodies["{2..12,28}(28=2*14)"]=list(range(2,13))+[28]
    bodies["{2..12,26}(26=2*13)+14?"]=None  # not covering, skip
    bodies["{2,3,4,6,7,8,9,12,11,13,14,10}"]=[2,3,4,6,7,8,9,10,11,12,13,14]
    bodies["even-heavy {2,4,6,8,10,12,14,3,9,5,11,13}"]=[2,4,6,8,10,12,14,3,9,5,11,13]
    bodies["{4,6,8,9,10,12,14,15,21,22,26,33}"]=[4,6,8,9,10,12,14,15,21,22,26,33]
    # add several {2..12,k} with k a small mult of 14 (backbone), and {3..14} type
    for k in [42,56,70,98,126,154,168,182,364]:
        bodies[f"{{2..12,{k}}}"]=list(range(2,13))+[k]
    worst=0.0; wb=None
    print(f"{'body B':<34} | cover({1}UB)? | coreCover | eps_1^2 = core*Sum eps^2 | vs 0.735")
    for name,B in bodies.items():
        if B is None: continue
        if len(B)!=12: continue
        full=[1]+B
        if not is_cov(full):
            print(f"{name:<34} | NOT covering, skip"); continue
        cc=coreCover(B); e2=(cc-1/7)**2
        if e2>worst: worst,wb=e2,name
        print(f"{name:<34} | {'yes':>10} | {cc:.4f} | {e2:.4f}                   | {'<0.735' if e2<tgt2 else 'OVER'}")
    print()
    import math
    print(f"WORST |core|=1: core*Sum eps^2 = {worst:.4f} at [{wb}]; sqrt = {math.sqrt(worst):.4f} (= Sum|eps| = coreCover-1/7)")
    print(f"MARGIN to opus's ceiling 0.735: {tgt2-worst:.4f} (square) / {6/7-math.sqrt(worst):.4f} (sqrt: 6/7={6/7:.4f})")
    print(f"=> opus's 0.328 is NOT the worst case (worst ~{worst:.2f}); the true binding |core|=1 margin to 6/7 is thin")
    print(f"   (~{6/7-math.sqrt(worst):.2f}). The tight large-sieve must clear {worst:.2f} < 0.735, not 0.328 -- less room than stated.")

if __name__ == "__main__":
    main()
