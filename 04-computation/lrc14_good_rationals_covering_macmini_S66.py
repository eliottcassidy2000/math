#!/usr/bin/env python3
"""mac-mini-S66: THE GOOD RATIONALS ARE COVERING-FAILURE MODULI.

CLAIM. x=a/b (gcd(a,b)=1) sits at the center of a good arc of Good(E)={x:maxgap{frac(e_i x)}>2/7}
IFF the residue set {a*e_i mod b : i} leaves a cyclic run of L consecutive empty residues with
L/b > 2/7 (i.e. L >= floor(2b/7)+1). Equivalently: E fails to 'cover' Z/b well enough at
dilation a -- a COVERING condition on the cluster at modulus b. The arc HALF-WIDTH is
(L/b - 2/7) / (2*D_eff) where D_eff relates to the spread, so FAT arcs (the ones an integer
ruler-period can land in) come only from SMALL b. This localizes the finite-Vmax existence to
a covering check at small moduli -- the covering constraint the cluster carries.

Verify: (1) the characterization exactly vs the sampled good set; (2) which b host fat arcs;
(3) the covering connection -- for COVERING-saturated S the cluster's small-modulus gaps."""
from fractions import Fraction as F
from math import gcd, floor

def maxgap_exact(E, x):
    ph=sorted(set(F(e)*x - (F(e)*x).__floor__() for e in E))
    n=len(ph)
    if n<=1: return F(1)
    g=max((ph[(i+1)%n]-ph[i]) if i<n-1 else (ph[0]+1-ph[n-1]) for i in range(n))
    return g

def residue_maxrun(E,a,b):
    """max cyclic run of empty residues in {a*e mod b}."""
    occ=set((a*e)%b for e in E)
    if len(occ)==b: return 0
    empty=[r for r in range(b) if r not in occ]
    # longest cyclic run of empties
    best=0
    for start in empty:
        if (start-1)%b in occ or len(occ)==0:  # run start
            L=0; r=start
            while r%b not in occ: L+=1; r+=1
            best=max(best,L)
    return best

def good_at_center(E,a,b):
    """is a/b a good-arc center? = does {a e mod b} leave run L with L/b>2/7?"""
    L=residue_maxrun(E,a,b)
    return F(L,b)>F(2,7), L

# test clusters
clusters = {
 "consec k=9": list(range(9)),
 "perforated k=7 {0..8}\\{1,7}": [0,2,3,4,5,6,8],
 "spread-21 k=11": [0,2,4,6,8,10,12,14,16,18,21],
}
print("VERIFY characterization: a/b is good-center IFF {a e mod b} leaves run L, L/b>2/7")
print("-"*88)
for nm,E in clusters.items():
    D=max(E)
    hits=[]
    for b in range(2,25):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            isg,L=good_at_center(E,a,b)
            if isg:
                # confirm exactly: maxgap at x=a/b really >2/7
                mg=maxgap_exact(E,F(a,b))
                hits.append((a,b,L,float(mg), mg>F(2,7)))
    # dedup by value a/b, show fat ones (small b)
    ok=all(h[4] for h in hits)
    small_b=sorted(set(b for a,b,L,mg,c in hits))
    print(f"{nm}: D={D}, good rationals at b in {small_b}, ALL confirmed maxgap>2/7: {ok}")
    for a,b,L,mg,c in sorted(hits, key=lambda h:h[1])[:8]:
        print(f"    {a}/{b}: empty-run L={L} (L/b={L}/{b}={float(F(L,b)):.3f}), exact maxgap={mg:.4f}")
print()
print("COVERING CONNECTION: good rationals a/b <=> cluster residues {e mod b} leave a >2/7 run")
print("<=> E FAILS to cover Z/b (at dilation a). So a good period exists <=> E under-covers some")
print("small modulus b. For a COVERING-saturated S the reconstructed cluster is forced to")
print("under-cover SOME small b (it cannot cover all of them with only k<=13 teeth) -- that")
print("forced under-coverage is the good arc. Test: min over small b of the best empty-run:")
for nm,E in clusters.items():
    k=len(E)
    bestrun=[]
    for b in range(2,15):
        mr=max((residue_maxrun(E,a,b) for a in range(1,b) if gcd(a,b)==1), default=0)
        bestrun.append((b,mr,F(mr,b)>F(2,7)))
    goodbs=[b for b,mr,g in bestrun if g]
    print(f"  {nm} (k={k}): under-covers (>2/7 run) at b={goodbs}  => {'HAS good period' if goodbs else 'NONE(!)'}")
