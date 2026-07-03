#!/usr/bin/env python3
"""
FAR-PEEL COMPLETENESS AUDIT (mac-mini-2026-07-03-S24): does opus's single-step peel + THM-609 close the HARD
covering families? Peel the LARGEST far runner w (>22), base B = other 12 (LRC(<=13) citation => good-region
floor THM-609). The family is lonely iff length(goodRegion(full)) > 0 = length(goodRegion B) - overlap(comb w).
Audit the hard cases: tight {1..12,182}, dilated APs, aligned band-blockers, GW-like, random covering.
A family FAILS the single-step peel iff removing the largest far comb kills the base good region.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))

def goodlen(speeds, Ng=600000):
    """measure of {t in [0,1): min_i ||v_i t|| >= 1/14}."""
    c=0; thr=1/14
    for k in range(Ng):
        t=(k+0.5)/Ng
        if min(nd(v*t) for v in speeds) >= thr: c+=1
    return c/Ng

def single_peel_ok(sp):
    """peel largest |v|>22; base = rest; family lonely iff full good region > 0. Return (baseGood, fullGood)."""
    A = sorted(sp, key=abs)
    w = A[-1]; base = A[:-1]
    return goodlen(base), goodlen(sp), w

if __name__ == "__main__":
    print("FAR-PEEL COMPLETENESS AUDIT: peel largest far runner; does base good region survive?")
    print("=" * 96)
    fams = {
        "tight {1..12,182}": list(range(1,13))+[182],
        "dilated AP x2 {2..26}": [2*i for i in range(1,14)],
        "dilated AP x3 (non-cov)": [3*i for i in range(1,14)],
        "aligned blocker ~1000": None,   # constructed below
        "aligned blocker ~30000": None,
        "GW-like doubling": [1,2,3,4,6,8,12,16,24,32,48,64,96],
        "random covering 1": None,
        "random covering 2": None,
    }
    rng = random.Random(24)
    def aligned(N):
        band=list(range(15,60)); rng.shuffle(band)
        far=sorted({q*round(N/q) for q in band[:10]}); far=[f for f in far if f>22]
        sp=far[:]
        for q in [8,9,5,7,11,13,2,3,4,6]:
            if len(sp)>=13: break
            if not any(s%q==0 for s in sp): sp.append(q)
        while len(sp)<13: sp.append(rng.randint(2,22))
        return sorted(set(sp))[:13]
    fams["aligned blocker ~1000"]=aligned(1000)
    fams["aligned blocker ~30000"]=aligned(30000)
    def randcov():
        for _ in range(2000):
            sp=sorted(rng.sample(range(1,400),13))
            if is_covering(sp) and any(v>22 for v in sp): return sp
        return None
    fams["random covering 1"]=randcov(); fams["random covering 2"]=randcov()

    print(f"{'family':>28} {'covering':>9} {'#far>22':>8} {'base good':>10} {'FULL good':>10} {'peel OK?':>9}")
    for name, sp in fams.items():
        if sp is None or len(set(sp))!=13:
            print(f"{name:>28} {'(skip)':>9}"); continue
        sp=sorted(set(sp)); cov=is_covering(sp); nf=sum(1 for v in sp if abs(v)>22)
        bg, fg, w = single_peel_ok(sp)
        ok = fg > 1e-6
        print(f"{name:>28} {str(cov):>9} {nf:>8} {bg:>10.5f} {fg:>10.6f} {str(ok):>9}")
    print("\n=> if FULL good > 0 for every COVERING family: the single-step far-peel (base LRC(13) floor THM-609")
    print("   minus the largest far comb) closes them -- opus's CoveringFarLonely 22 route is complete on these.")
    print("   (non-covering families are lonely at t=1/q, q<=14, separately -- not the peel's job.)")
