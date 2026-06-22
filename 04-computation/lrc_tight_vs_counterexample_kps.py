"""
CRITICAL CHECK (kps-S31n): for each safe-measure-0 set, is it TIGHT (M=1/14, witness exists)
or a COUNTEREXAMPLE (M<1/14, safe set empty)?  M = max_t min_s ||s t||.
The max-min is achieved at a rational t = a/(L) where two ||s t|| balance; candidate t are the
breakpoints (14m+-1)/(14 s).  Evaluate min_s ||s t|| there and over a fine grid; report M.
"""
from fractions import Fraction as F

def nfrac(x):  # ||x|| as Fraction in [0,1/2]
    r = x % 1
    return min(r, 1-r)

def lonely_constant(S, grid=200000):
    S=[s for s in S if s!=0]
    # candidate witness points: breakpoints (14m+-1)/(14 s) AND midpoints, plus a fine grid
    cands=set()
    for s in S:
        a=abs(s)
        for m in range(0, a+1):
            for sgn in (-1,1):
                t=F(14*m+sgn, 14*a)
                if 0<=t<=1: cands.add(t)
    # also j/(n+1) style
    for d in range(1, 60):
        for j in range(0, d+1):
            cands.add(F(j,d))
    best=F(0); arg=None
    for t in cands:
        mn=min(nfrac(s*t) for s in S)
        if mn>best: best=mn; arg=t
    # fine grid double check (float)
    bestf=0.0; argf=0.0
    for i in range(1,grid):
        t=i/grid
        mn=min(min((s*t)%1, 1-(s*t)%1) for s in S)
        if mn>bestf: bestf=mn; argf=t
    return best, arg, bestf, argf

for name,S in [("AP {1..13}", list(range(1,14))),
               ("GW {1..11,13,24}", list(range(1,12))+[13,24]),
               ("2*AP", [2*x for x in range(1,14)])]:
    M, arg, Mf, argf = lonely_constant(S)
    status = "TIGHT (witness exists, LRC holds)" if M>=F(1,14) else "*** COUNTEREXAMPLE M<1/14 ***"
    print(f"{name:20s}: M_exact={M}={float(M):.6f} at t={arg};  grid M={Mf:.6f} at {argf:.5f}")
    print(f"  {'':18s}  vs 1/14={1/14:.6f}  =>  {status}")
