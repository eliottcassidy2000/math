#!/usr/bin/env python3
"""
lonely_set_equidistribution_klein.py  --  klein-2026-07-01-S65

EQUIDISTRIBUTION ON THE FIXED LONELY SET L_C. For the covering-min construction C={1..n-2,n(n-1)}, the
lonely set at level r is  L_C(r) = { t in [0,1) : min_{v in C} ||v t|| > r }  (times where the observer
is lonely).  A BEATER must reduce M below M_C by having speeds whose danger arcs cover L_C(r). A large
speed w helps iff its danger {t: ||w t|| < r'} COVERS a lot of L_C.

KEY QUESTION (extending S64's far-element/equidistribution finding): does {w t mod 1} EQUIDISTRIBUTE on
L_C (so w's danger covers only ~2r'*meas(L_C) of it -- inefficient, can't cheaply kill L_C), or can a
TUNED / RESONANT large w CONCENTRATE its danger on L_C (the CRT-escape)? By Weyl, generic large w
equidistributes; resonance (arithmetic structure of L_C) is the escape.

Computes: (1) meas + interval-count of L_C(r) (its Cantor/three-gap structure); (2) the coverage
frac_w = meas(L_C ∩ {||wt||<r'}) / meas(L_C) for w=1..W; the equidistribution prediction 2r'; the
histogram (concentration at 2r' = equidistribution; outliers = resonances = CRT-escape candidates).
"""
import numpy as np
from fractions import Fraction as F

def norms(v, t):  # ||v t|| vectorized, t in [0,1)
    x = (v*t) % 1.0
    return np.minimum(x, 1-x)

if __name__=="__main__":
    n=14; C=list(range(1,n-1))+[n*(n-1)]; Mc=F(n,n*n-n+1)  # 14/183
    N=400000; t=np.arange(N)/N
    # min over construction speeds
    G = np.full(N, 1.0)
    for v in C: G = np.minimum(G, norms(v,t))
    print(f"construction C={C}, M_C = {Mc} = {float(Mc):.5f}")
    print("(1) THE FIXED LONELY SET L_C(r) = {t: min_v||vt|| > r}: measure + #intervals")
    for r in [0.05, 0.06, 0.07, 1/n, 13/183, float(Mc)-1e-6]:
        mask = G > r
        meas = mask.mean()
        # count intervals (runs of True)
        d = np.diff(mask.astype(int)); nint = int((d==1).sum()) + (1 if mask[0] else 0)
        print(f"   r={r:.5f}: meas(L_C)={meas:.5f}  #intervals={nint}  (r vs floor 1/n={1/n:.5f}, M_C={float(Mc):.5f})")

    # fix a level for the equidistribution test: just below M_C (the near-binding lonely set)
    r = float(Mc) - 0.002
    Lmask = G > r; Lmeas = Lmask.mean()
    rp = float(Mc)   # a new speed's danger radius (to reach the covmin level)
    print()
    print(f"(2) EQUIDISTRIBUTION on L_C(r={r:.4f}) (meas={Lmeas:.5f}): coverage frac_w of L_C by w's danger")
    print(f"    {{t: ||wt||<r'}}, r'={rp:.4f}; Weyl prediction 2r'={2*rp:.4f}. Concentration=equidistrib; outliers=resonance.")
    W=1200; fracs=[]
    for w in range(1, W+1):
        cov = ( (norms(w,t) < rp) & Lmask ).mean() / Lmeas
        fracs.append((cov, w))
    fr = np.array([f for f,_ in fracs])
    pred = 2*rp
    print(f"    mean frac_w = {fr.mean():.4f} (vs 2r'={pred:.4f}); std={fr.std():.4f}")
    top = sorted(fracs, reverse=True)[:8]; bot = sorted(fracs)[:8]
    print(f"    MOST-covering w (resonant, CRT-escape candidates): {[(round(f,3),w) for f,w in top]}")
    print(f"    LEAST-covering w (anti-resonant, lonely-preserving): {[(round(f,3),w) for f,w in bot]}")
    # is the construction's own killer 182 resonant or anti? and multiples of Phi6?
    for w in [182, 183, 13, 14, 7, 61, 183*2]:
        cov=((norms(w,t)<rp)&Lmask).mean()/Lmeas
        print(f"    w={w:>4}: frac={cov:.4f}  ({'RESONANT (covers L_C)' if cov>pred*1.3 else 'ANTI (misses L_C)' if cov<pred*0.7 else 'equidistributed'})")
    print()
    print("=> if generic large w -> frac ~ 2r' (equidistribution) but a FEW resonant w cover much more,")
    print("   the CRT-escape = w RESONATING with the arithmetic (three-gap/Farey) structure of L_C.")
    print("   Extending S64: the far element beats ONLY by resonating with L_C; generic large speeds")
    print("   equidistribute and waste their danger off the lonely set.")
