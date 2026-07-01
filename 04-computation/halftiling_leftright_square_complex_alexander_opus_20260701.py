"""Is the half-tiling a LEFT-RIGHT square complex? The tiling cube (Z/2)^m carries Z/2xZ/2 = <flip, sigma>:
flip = complement-tiling (LEFT), sigma = grid-reflection TRANS (RIGHT), commuting => SQUARES {t,flip t,sigma t,
flip sigma t}, degenerating to LINES (size 2) at sigma-fixed (grid-sym=blue=half-tiling) & flipsigma-fixed.
Plus Alexander duality b0(lonely)=b0(danger) on S^1 (kps-S22)."""
import numpy as np
from itertools import product
def sq_structure(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES); ti={t:i for i,t in enumerate(TILES)}
    TRANS=[ti[(n-y+1,n-x+1)] for (x,y) in TILES]; full=(1<<m)-1
    def flip(t): return t^full
    def sig(t):
        r=0
        for i in range(m):
            if (t>>i)&1: r|=1<<TRANS[i]
        return r
    orbits=[]; seen=set()
    for t in range(1<<m):
        if t in seen: continue
        orb={t,flip(t),sig(t),flip(sig(t))}
        for x in orb: seen.add(x)
        orbits.append(len(orb))
    from collections import Counter
    c=Counter(orbits)
    gridsym=sum(1 for t in range(1<<m) if sig(t)==t)   # sigma-fixed = half-tiling = 2^D
    flipsig=sum(1 for t in range(1<<m) if flip(sig(t))==t)  # sigma=complement-tiling
    return m,dict(c),gridsym,flipsig
print("Z/2xZ/2 = <flip,sigma> square structure on the tiling cube:")
for n in [4,5,6]:
    m,c,gs,fs=sq_structure(n)
    print(f"  n={n} (m={m}, 2^m={1<<m} tilings): orbit sizes {c}  (size4=SQUARES, size2=LINES); sigma-fixed(half-tiling/blue)={gs}=2^{gs.bit_length()-1}, flipsigma-fixed={fs}")
    sq=c.get(4,0); ln=c.get(2,0); print(f"     => {sq} squares + {ln} lines; check 4*{sq}+2*{ln}={4*sq+2*ln}=2^m={1<<m}: {4*sq+2*ln==(1<<m)}")
# Alexander duality on S^1: b0(lonely)=b0(danger)
def bettis(S,r,N=600000):
    t=np.arange(N)/N; keep=np.ones(N,bool)
    for v in S:
        keep &= (np.abs(((v*t+0.5)%1)-0.5) >= r-1e-12)
    def comps(mask):
        prev=mask[-1]; c=0
        for x in mask:
            if x and not prev: c+=1
            prev=x
        return c
    return comps(keep), comps(~keep)
print("\nAlexander duality on S^1 (lonely = complement of danger): b0(lonely) == b0(danger)?")
for n in [5,7,14]:
    S=list(range(1,n)); bl,bd=bettis(S,0.99/n)
    print(f"  AP n={n} (r=0.99/n): b0(lonely)={bl}, b0(danger)={bd}  => equal? {bl==bd} (Alexander dual on S^1)")
