#!/usr/bin/env python3
"""
Is the measure-zero good-region GAP only gcd>1 dilations? (mac-mini-2026-07-03-S24)
opus's CoveringFarLonely (no gcd=1) includes {2,4,..,26} (covering, far entry, but good region = 0: lonely
only at isolated t=1/28, M=1/14 EXACTLY). The positive-length peel (THM-609 + far_peel) needs length>0, so it
CANNOT reach measure-zero families. FIX: gcd-reduce (dilation invariance). This holds iff EVERY covering
far-entry family with measure-zero good region is a gcd>1 dilation of a window/non-covering family, i.e. NO
gcd=1 covering far-entry family is measure-zero. Search for a gcd=1 counterexample.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def gcd_all(xs): return reduce(gcd, xs)
def nd(x): x = x % 1; return min(x, 1-x)
def is_covering(sp): return all(any(v % q == 0 for v in sp) for q in range(2,15))

def M_and_goodlen(speeds, Nscan=400000, Ng=400000):
    """M = max_t min_i ||v_i t|| (fine scan); goodlen = measure{min>=1/14}."""
    vmax=max(speeds); K=25; s=min(Nscan, vmax*K); M=0.0
    for k in range(1,s):
        t=k/s; m=min(nd(v*t) for v in speeds)
        if m>M: M=m
    thr=1/14; c=0
    for k in range(Ng):
        t=(k+0.5)/Ng
        if min(nd(v*t) for v in speeds)>=thr: c+=1
    return M, c/Ng

if __name__ == "__main__":
    rng = random.Random(77)
    print("Searching for a gcd=1 COVERING far-entry family with MEASURE-ZERO good region (M ~ 1/14).")
    print("=" * 88)
    print(f"1/14 = {1/14:.6f}; 14/183 = {float(F(14,183)):.6f} (tight covering-min, positive good region)")
    print("-" * 88)
    # (A) confirm the gap examples
    for name, sp in [("{2..26}=2*{1..13}", [2*i for i in range(1,14)]),
                     ("{3..39}=3*{1..13}", [3*i for i in range(1,14)]),
                     ("{1..12,182} tight g=1", list(range(1,13))+[182])]:
        M,g = M_and_goodlen(sp)
        print(f"  {name:>24}: gcd={gcd_all(sp)}, covering={is_covering(sp)}, far={sum(1 for v in sp if v>22)}, M={M:.6f}, goodlen={g:.6f}")
    print("-" * 88)
    # (B) hunt gcd=1 covering far-entry families with tiny good region (near-tight)
    worst = (1.0, None)   # smallest goodlen found
    n_tested = 0; n_zero = 0
    for _ in range(3000):
        # bias toward near-tight: APs, near-APs, {1..k}+outliers
        kind = rng.choice(["ap","window+outlier","near-tight","random"])
        if kind=="ap":
            a=rng.randint(1,3); d=rng.randint(1,3)
            sp=[a+d*i for i in range(13)]
        elif kind=="window+outlier":
            base=rng.sample(range(1,15), rng.randint(9,12))
            sp=sorted(set(base))+[rng.randint(23,400) for _ in range(13-len(set(base)))]
            sp=sorted(set(sp))
        elif kind=="near-tight":
            sp=list(range(1,13))+[rng.choice([182,183,181,89,155,167])]
        else:
            sp=sorted(rng.sample(range(1,300),13))
        if len(set(sp))!=13: continue
        sp=sorted(set(sp))
        if gcd_all(sp)!=1 or not is_covering(sp) or not any(v>22 for v in sp): continue
        n_tested+=1
        M,g = M_and_goodlen(sp, Nscan=150000, Ng=150000)
        if g < worst[0]: worst=(g, sp, M)
        if g < 1e-5: n_zero+=1
    print(f"tested {n_tested} gcd=1 covering far-entry families")
    print(f"SMALLEST good-region length found: {worst[0]:.7f}  (family {worst[1]}, M={worst[2]:.6f})")
    print(f"# with goodlen < 1e-5 (measure-zero): {n_zero}")
    print(f"\n=> if smallest goodlen > 0 (no measure-zero): gcd=1 covering far families ALL have positive good")
    print("   region => THM-609 + far-peel reaches them; the measure-zero gap is ONLY gcd>1 dilations,")
    print("   removed by WLOG gcd=1 (dilation invariance). opus's CoveringFarLonely NEEDS the gcd=1 hypothesis")
    print("   (or a dilation-reduction lemma) to be provable by the positive-length route.")
