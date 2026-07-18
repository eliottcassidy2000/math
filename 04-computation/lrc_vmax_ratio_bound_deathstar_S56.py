from fractions import Fraction as F
from math import gcd
n=14

def good_components(fam):
    # components (as (lo,hi) fractions) of {t in [0,1): min_v ||v t|| > 1/n}, exact.
    # danger of v: intervals [(n k -1)/(n v), (n k +1)/(n v)] for k=0..v (wrap). Merge, complement.
    danger=[]
    for v in fam:
        for k in range(0,v+1):
            lo=F(n*k-1,n*v); hi=F(n*k+1,n*v)
            danger.append((lo,hi))
    # clip to [0,1) with wrap: just work on [0,1] and also add wrapped copies
    ivs=[]
    for lo,hi in danger:
        ivs.append((lo,hi))
    ivs.sort()
    # merge
    merged=[]
    for lo,hi in ivs:
        if merged and lo<=merged[-1][1]:
            merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else:
            merged.append((lo,hi))
    # good components = gaps in [0,1] between merged danger (danger covers near 0 and 1 by k=0,v)
    comps=[]
    prev=F(0)
    for lo,hi in merged:
        if lo>prev:
            comps.append((prev,lo))
        prev=max(prev,hi)
    if prev<F(1):
        comps.append((prev,F(1)))
    # widths
    return [(a,b,b-a) for a,b in comps if b>a]

def Vmax_check(name,fam):
    fam=sorted(fam); w=fam[-1]; v2=fam[-2]; Vp=fam[:-1]
    comps=good_components(Vp)
    Lmax=max((c[2] for c in comps), default=F(0))
    cov_bound=F(2,n*w)           # each component must be <= this
    lower=F(2,n*(n-1)*v2)        # provable lower bound on largest component
    print(f"{name}: Vmax=w={w}, v2={v2}")
    print(f"  #good-components of V\\{{w}} = {len(comps)}; Lmax={Lmax}={float(Lmax):.6f}")
    print(f"  covering bound 2/(n w)={cov_bound}={float(cov_bound):.6f}; all comps <= it? {all(c[2]<=cov_bound for c in comps)}")
    print(f"  provable lower 2/(n(n-1)v2)={float(lower):.6f} <= Lmax? {lower<=Lmax}")
    print(f"  => BOUND  Vmax <= (n-1)*v2 = {(n-1)*v2}:  {w} <= {(n-1)*v2}  {'OK' if w<=(n-1)*v2 else 'VIOLATED(=>not tight)'}")
    print()

Vmax_check("AP {1..13} (tight)", list(range(1,14)))
Vmax_check("GW {1..11,13,24} (tight)", [1,2,3,4,5,6,7,8,9,10,11,13,24])
Vmax_check("deep well {1..12,182} (NON-tight, M=14/183)", list(range(1,13))+[182])
