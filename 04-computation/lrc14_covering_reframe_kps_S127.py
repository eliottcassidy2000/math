# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont64: work the covering reframe -- can 12 smooth 1/14-collars tile [1/14, 13/14]?
#
# cont.63 reframe: {1}uB lonely (LRC) <=> the 12 non-core collars C_i={t:||b_i t||<1/14} leave a GAP in runner-1's
# good region [1/14,13/14]. This tests whether the reframe gives an OBSTRUCTION (measure/structure) or is LRC-hard.
from math import gcd
from fractions import Fraction as F
def norm(x): r=x-int(x); r=r+1 if r<0 else r; return min(r,1-r)
def is_cov(v,N=14): return all(any(x%d==0 for x in v) for d in range(2,N+1))

def collar_measure_in_interval(B, level=1/14, a=1/14, b=13/14, grid=800000):
    """measure of (union of collars) INSIDE [a,b], and total collar measure, and the GAP measure (uncovered)."""
    cov=0; gap=0
    for k in range(grid):
        t=a+(b-a)*k/grid
        bad=any(norm(bi*t)<level for bi in B)
        if bad: cov+=1
        else: gap+=1
    w=(b-a)
    return cov/grid*w, gap/grid*w

def near_endpoint_covering(B, level=1/14, eps_grid=2000):
    """which runners cover just above 1/14? claim: only mult-of-14 (arc at 1/14) reach there."""
    reaches=[]
    for bi in B:
        # bi covers t just above 1/14 iff its nearest arc-center j/bi is within level/bi of 1/14
        # i.e. |bi/14 - round(bi/14)| < level  (arc at 1/14 exists iff bi*1/14 near integer)
        r = norm(bi*(1/14))
        if r < level:  # bi is bad AT 1/14 -> covers a nbhd; reach above = 1/14 + (level - r)/bi
            reaches.append((bi, 1/14 + (level - r)/bi))
    return reaches

def main():
    print("THE COVERING REFRAME: {1}uB lonely <=> 12 collars leave a GAP in [1/14,13/14].")
    print(f"interval width = 12/14 = {12/14:.5f}; each collar measure = 2/14 = 1/7; total 12 collars = 12/7 = {12/7:.4f}\n")
    print("(1) MEASURE budget -- does the interval fit the collars? (covered vs gap measure inside [1/14,13/14]):")
    bodies = {
        "deep well {2..12,182}": list(range(2,13))+[182],
        "ladder {2..12,364}": list(range(2,13))+[364],
        "multi {2..11,13,84,?}": [2,3,4,5,6,7,8,9,10,11,13,84],
        "{3..14}": list(range(3,15)),
    }
    for name,B in bodies.items():
        cov,gap = collar_measure_in_interval(B)
        print(f"    {name:<24} covered={cov:.5f}  GAP={gap:.5f}  (gap>0 => LRC holds; gap is the poke-out)")
    print("    => total collar measure 12/7=1.71 >> interval 0.857: MEASURE PERMITS covering. No measure obstruction.")

    print("\n(2) ENDPOINT structure -- just above 1/14, ONLY mult-of-14 runners cover (arc centred at 1/14):")
    for name,B in [("deep well body", list(range(2,13))+[182]), ("ladder S_2 body", list(range(2,13))+[364])]:
        r = near_endpoint_covering(B)
        print(f"    {name:<16}: runners bad at 1/14 (cover a nbhd) = {[(b,round(reach,5)) for b,reach in r]}")
        if r:
            m14 = [b for b,_ in r]
            print(f"                     (only mult-of-14 = {m14}; each covers up to 1/14 + (1/14)/b = tiny reach)")

    print("\n(3) the GAP LOCATION varies (no universal localization): single-killer near 1/14, multi-killer mid-interval:")
    from math import gcd as g
    def Mexact(v):
        qc=4*max(v)+2; best=F(0); bt=None
        for q in range(2,qc):
            for p in range(1,q):
                if g(p,q)==1:
                    m=min(norm(F(vi*p,q)) for vi in v)
                    if m>best: best,bt=m,F(p,q)
        return best,bt
    for name,v in [("single-killer {1..12,182}", list(range(1,13))+[182]),
                   ("multi-killer {1..11,13,84}", list(range(1,12))+[13,84])]:
        M,t = Mexact(v)
        print(f"    {name:<26} gap (lonely t*) = {t} = {float(t):.4f}  (||t*||={float(norm(F(t.numerator,t.denominator))):.4f})")
    print()
    print("=> THE REFRAME IS LRC-HARD: measure permits covering (no budget obstruction); the gap location varies")
    print("   (near-1/14 single-killer, mid-interval multi-killer, no universal gap point); the endpoint structure")
    print("   is only a tiny-nbhd fact. So the covering reframe is EQUIVALENT to the |core|=1 residual, NOT a")
    print("   shortcut -- the crux stays the fine loneliness/discrepancy (opus Fourier). The single-killer slice")
    print("   IS proved (in covering language: {2..12,182c} collars leave the gap 14c/(182c+1), Lean cont.60/61).")

if __name__ == "__main__":
    main()
