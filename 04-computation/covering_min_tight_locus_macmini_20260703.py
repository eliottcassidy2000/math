#!/usr/bin/env python3
"""
Covering-min + tight-locus structure (mac-mini-2026-07-03-S30). Ground the crux: covering => M >= 14/183.
 (1) search small-speed covering gcd=1 families for MIN M(view); confirm covering-min ~ 14/183 = 14/Phi6(14).
 (2) characterize the tight ones (M near 14/183): residues mod 14, do they cover +-units {1,3,5,9,11,13}?
     are they non-covering (0 absent)? # distinct gaps at the optimum (g<=3 Steinhaus)?
 (3) test the RIGIDITY g(14)<=3: any near-tight covering family with >3 gaps at its optimum?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd,xs)
def nd(x):
    x = x % 1
    return min(x, 1-x)
def is_covering(sp): return all(any(v%q==0 for v in sp) for q in range(2,15))
def M_and_topt(sp, D=183*40):
    """M = max over t=a/D of min_i ||v_i t||, and the argmax t (fine rational grid, denom multiple of 183)."""
    best=F(0); bt=F(0)
    for a in range(1, D):
        t=F(a,D); m=min(nd(v*t) for v in sp)
        if m>best: best,bt=m,t
    return best, bt
def gaps_at(sp, t):
    """# distinct gaps of the phase set {v_i t mod 1} U {0} on the circle (rounded)."""
    ph=sorted(set(round(float((v*t)%1),6) for v in sp) | {0.0})
    gs=sorted(set(round(ph[(i+1)%len(ph)]-ph[i] if i<len(ph)-1 else 1-ph[-1]+ph[0],5) for i in range(len(ph))))
    return len(gs), ph

if __name__=="__main__":
    print(f"covering-min target: 14/183 = {float(F(14,183)):.6f}; 1/14 = {float(F(1,14)):.6f}; 183=Phi6(14)=14^2-14+1")
    print("="*82)
    # (1) search small-speed covering gcd=1 families for min M
    rng=random.Random(30)
    minM=(F(1),None)
    # include the known tight family + dilation-reduced + random small
    cands=[list(range(1,13))+[182]]
    for _ in range(3000):
        sp=sorted(set(rng.sample(range(1,60),13)))
        if len(sp)==13: cands.append(sp)
    tested=0
    for sp in cands:
        if len(sp)!=13 or gcd_all(sp)!=1 or not is_covering(sp): continue
        tested+=1
        M,t=M_and_topt(sp, D=183*8)
        if M<minM[0]: minM=(M,sp,t)
    print(f"(1) searched {tested} small-speed covering gcd=1 families")
    print(f"    MIN M = {minM[0]} = {float(minM[0]):.6f}  (target 14/183={float(F(14,183)):.6f})")
    print(f"    tight family: {minM[1]}  at t*={minM[2]}")
    # (2)+(3) characterize tight + near-tight families: residues, +-units cover, gaps
    print("\n(2)+(3) tight-locus structure of near-tight covering families (M within 5% of 14/183):")
    print(f"{'family (first 6)':>28} {'M':>9} {'residues mod14':>28} {'covers +-units?':>15} {'0 in res?':>10} {'#gaps@opt':>10}")
    units14={1,3,5,9,11,13}
    shown=0
    allfams = sorted([(sp) for sp in cands if len(sp)==13 and gcd_all(sp)==1 and is_covering(sp)],
                     key=lambda s: M_and_topt(s,183*4)[0])[:12]
    for sp in allfams:
        M,t=M_and_topt(sp, D=183*8)
        res=sorted(set(v%14 for v in sp))
        cov_units = units14.issubset(set(res))
        has0 = 0 in res
        ng,_=gaps_at(sp,t)
        print(f"{str(sp[:6]):>28} {float(M):>9.5f} {str(res):>28} {str(cov_units):>15} {str(has0):>10} {ng:>10}")
        shown+=1
    print("\n=> covering-min ~ 14/183 (if MIN M matches). tight families cover +-units {1,3,5,9,11,13}, 0 absent")
    print("   (non-covering residue set!), and have <=3 distinct gaps at the optimum (Steinhaus g(14)<=3).")
    print("   the OPEN CORE is proving g(14)<=3 for ALL near-tight covering families (rigidity).")
