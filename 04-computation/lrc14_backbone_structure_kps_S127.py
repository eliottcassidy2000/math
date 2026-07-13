# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont65: push the endpoint mult-of-14 structure -> the BACKBONE + FILLERS decomposition.
#
# cont.64 endpoint fact: near 1/14 only the mult-of-14 runner covers. GLOBALIZE it: the mult-of-14 runner 14m
# is a BACKBONE -- it is bad at EVERY grid point j/14 (14m*(j/14)=jm in Z) AND has intermediate arcs at
# k/(14m), tiling [1/14,13/14] on a grid of spacing 1/(14m), leaving inter-arc GAPS of width ~13/(196m). The
# other 11 runners are FILLERS. The lonely point t* (the covering-min gap) must lie in a backbone gap (14m is
# lonely there). This tests: (1) t* is always in a backbone gap; (2) which fillers BIND at t*; (3) the
# grid-point modular covering law (near j/14, covered by mult-of-(14/gcd(j,14))).
from math import gcd
from fractions import Fraction as F
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def normF(x):
    x=x-(x.numerator//x.denominator);
    if x<0: x+=1
    return min(x,1-x)
def Mexact(v):
    qc=4*max(v)+2; best=F(0); bt=None
    for q in range(2,qc):
        for p in range(1,q):
            if gcd(p,q)==1:
                m=min(normF(F(vi*p,q)) for vi in v)
                if m>best: best,bt=m,F(p,q)
    return best,bt
def mult14(v): return [b for b in v if b%14==0]

def main():
    print("(1)+(2) the lonely point t* sits in a mult-of-14 BACKBONE gap; which runners BIND there:")
    print(f"{'family':<26} | M | t* | backbone (mult14) | binding runners at t* (dist=M)")
    fams = {
        "single-killer {1..12,182}": list(range(1,13))+[182],
        "ladder S_2 {1..12,364}": list(range(1,13))+[364],
        "multi-killer {1..11,13,84}": list(range(1,12))+[13,84],
        "multi {1..10,13,22,84}": list(range(1,11))+[13,22,84],
    }
    for name,v in fams.items():
        M,t = Mexact(v)
        bb = mult14(v)
        # binding runners: those with ||b t*|| == M
        binders = [b for b in v if normF(F(b*t.numerator, t.denominator))==M]
        # is 14m lonely at t* (in its gap)? yes iff ||14m t*|| >= M i.e. not covering
        bb_lonely = all(normF(F(b*t.numerator,t.denominator))>=M for b in bb)
        print(f"{name:<26} | {str(M):>7} | {str(t):>7} | {bb} | {binders}   (backbone lonely@t*: {bb_lonely})")
    print()
    print("=> single-killer: ONLY runner 1 + the backbone bind (=> the 2-runner balance, PROVED cont.60/61).")
    print("   multi-killer: MORE runners bind (the fillers engage) => the backbone gap is contested, mid-interval.")

    print("\n(3) the GRID-POINT modular covering law: near j/14, a runner b covers iff 14|bj iff b==0 mod 14/gcd(j,14):")
    print(f"    {'j':>2} | gcd(j,14) | need b == 0 mod | covered by")
    for j in range(1,14):
        g=gcd(j,14); need=14//g
        who = {2:"even (mult 2)",7:"mult 7",14:"mult 14"}.get(need,f"mult {need}")
        print(f"    {j:>2} | {g:>9} | {need:>14} | {who}")
    print("    => all 13 grid points are COVERED by any covering family (has mult 14 -> covers the 6 coprime-14")
    print("       points and (14=2*7) the 6 even + the j=7 point). So gaps are strictly BETWEEN grid points,")
    print("       inside the backbone's inter-arc slots -- the covering-min lives in one backbone gap.")

if __name__ == "__main__":
    main()
