#!/usr/bin/env python3
"""
The covering-min is a CONTINUED-FRACTION / Ostrowski object; the extremizer's gaps = classical THREE-GAP.
(mac-mini-2026-07-04-S38) Creative Fibonacci/Zeckendorf connection to the remaining core (covering-min/g<=3).
 (1) deep well {1..12,182} at t*=14/183: phases = {k*alpha mod 1}? gaps <=3? alpha=14/183=[0;13,14]?
 (2) the Ostrowski ladder M_k=[0;n-1,k]=k/(k(n-1)+1); covering-min = k=n (deep well). Zeckendorf=golden [0;1,1,..].
 (3) THREE-GAP THEOREM check: for alpha=n/Phi6, {k alpha mod 1: k=1..n-1} has <=3 gaps (Steinhaus) => g(n)<=3
     for the extremizer FROM the classical theorem (via the CF structure). Does it extend beyond the extremizer?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce

def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def cf(fr):
    """continued fraction of a Fraction in [0,1)."""
    a = []; x = fr
    for _ in range(12):
        if x == 0: break
        r = 1/x; ai = r.numerator//r.denominator; a.append(ai); x = r - ai
    return a
def gaps_of(points):
    """distinct gap lengths of a set of points (Fractions in [0,1)) on the circle."""
    p = sorted(set(points))
    g = set()
    for i in range(len(p)):
        nxt = p[(i+1) % len(p)]
        d = (nxt - p[i]) % 1 if i < len(p)-1 else (1 - p[-1] + p[0])
        g.add(d)
    # recompute cleanly (circular)
    g = set()
    for i in range(len(p)-1): g.add(p[i+1]-p[i])
    g.add(1 - p[-1] + p[0])
    return sorted(g)

if __name__ == "__main__":
    n = 14; phi6 = n*n-n+1  # 183
    print(f"n={n}, Phi6(n)={phi6}, covering-min = n/Phi6 = {F(n,phi6)} = {float(F(n,phi6)):.6f}")
    print(f"  continued fraction of n/Phi6: {cf(F(n,phi6))}   (expect [0; n-1, n] = [0;{n-1},{n}])")
    print()
    # (1) deep well at t*=14/183
    DW = list(range(1,13))+[182]
    t = F(n, phi6)
    phases = sorted(set((v*t) % 1 for v in DW) | {F(0)})
    print(f"(1) deep well {{1..12,182}} at t*={t}: phases (numerators over {phi6}):")
    print("   ", [ (v*t % 1).numerator*phi6//(v*t%1).denominator if (v*t%1)!=0 else 0 for v in sorted(DW)])
    g = gaps_of(phases)
    print(f"   phase set has {len(g)} distinct gaps: {[str(x)+f'({x.numerator})' for x in g]}  => g<=3? {len(g)<=3}")
    # compare to {k*alpha mod 1}: multiples of n mod Phi6
    core = sorted(set(F((k*n) % phi6, phi6) for k in range(1,n)) | {F(0)})
    gc = gaps_of(core)
    print(f"   {{k*n mod Phi6 : k=1..{n-1}}} (={{k*alpha}}, alpha=n/Phi6): {len(gc)} gaps {[x.numerator for x in gc]} => classical THREE-GAP")
    print()
    # (2) Ostrowski ladder
    print("(2) Ostrowski ladder M_k=[0;n-1,k]=k/(k(n-1)+1):")
    for k in [1,2,7,13,14,20]:
        Mk = F(k, k*(n-1)+1)
        print(f"   k={k}: M_k = {Mk} = {float(Mk):.5f}   (k={n} = covering-min {F(n,phi6)})")
    print(f"   golden/Zeckendorf limit: [0;1,1,1,...] = (sqrt5-1)/2 = 0.618 (all-1s CF); LRC covmin = 2-term [0;n-1,n].")
    print()
    # (3) three-gap for the extremizer vs a generic covering family
    print("(3) THREE-GAP at the optimum: extremizer (AP-like, {k alpha}) vs generic covering family:")
    fams = [("deep well {1..12,182}", DW), ("AP {1..13}", list(range(1,14))),
            ("generic cov {2,3,4,6,7,8,9,10,11,12,13,14,84}", [2,3,4,6,7,8,9,10,11,12,13,14,84])]
    def M_topt(sp, D=phi6*6):
        best=F(0); bt=None
        for a in range(1,D):
            tt=F(a,D); m=min(nd(v*tt) for v in sp)
            if m>best: best,bt=m,tt
        return best,bt
    for name, S in fams:
        M,tt = M_topt(S)
        ph = sorted(set((v*tt)%1 for v in S)|{F(0)})
        g = gaps_of(ph)
        # are the phases {k*alpha}? check if S's residues at tt form an AP
        print(f"   {name[:30]:>30}: M={float(M):.5f} t*={tt}  #gaps@opt={len(g)}  <=3? {len(g)<=3}")
    print("\n=> covering-min = [0;n-1,n] (Ostrowski/CF); extremizer phases = {k*alpha} => classical THREE-GAP g<=3.")
    print("   The g(n)<=3 rigidity FOLLOWS from Steinhaus IF tight => {k alpha}-config (=AP-like) = GAP-A itself.")
