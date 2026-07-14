#!/usr/bin/env python3
"""mac-mini-S89b: does the conc forbidden band STAY OPEN as Vmax grows? (uniform gap = last bit closed
on this stratum; band closing = need finer analysis). Track sup_covering conc vs Vmax, identify the
edge family, and contrast bounded-near-AP vs far (deep-well) families."""
import numpy as np
from itertools import combinations
c=1.0/14; N=400_000; t=(np.arange(N)+0.5)/N; band=(t<1.0/14)
def conc(C):
    ok=np.ones(N,bool)
    for w in C:
        r=(w*t)%1.0; ok &= (np.minimum(r,1-r)>=c)
    mg=ok.mean()
    return (14*(ok&band).mean()/mg, mg) if mg>0 else (None,0.0)
def covers(S): return all(any(v%q==0 for v in S) for q in range(2,15))
print("sup covering conc vs Vmax (edge family = closest covering set to the AP wall 7):")
for Vmax in [15,16,17,18]:
    pool=list(range(2,Vmax+1)); sup=-1; arg=None
    for C in combinations(pool,12):
        S=(1,)+C
        if max(C)!=Vmax: continue      # only NEW families at this Vmax
        if not covers(S): continue
        a,mg=conc(list(C))
        if a is not None and a>sup: sup=a; arg=S
    if arg: print(f"  Vmax={Vmax}: sup conc={sup:.4f}  margin(1-sup/7)={1-sup/7:.4f}  edge={arg}")
# explicit: drop-x near-AP family conc (the band edge), all x
print("\nnear-AP bounded covering {1..14}\\{x} -- conc by dropped speed x (the apex structure):")
for x in range(2,8):
    S=tuple(v for v in range(1,15) if v!=x); C=[v for v in S if v!=1]
    a,mg=conc(C); print(f"  drop-{x}: conc={a:.4f}  L>={mg*(1-a/7):.5f}   ({'closest to AP wall' if x==6 else ''})")
# contrast: far families (deep well, near-AP residue) -- conc much lower (further from AP)
print("\nfar families (handled by disc_v/THM-733-5) -- conc far below the band edge:")
for name,S in [("deep well {1..12,182}",tuple(list(range(1,13))+[182])),
               ("near-AP resid {1..11,13,84}",tuple(list(range(1,12))+[13,84]))]:
    C=[v for v in S if v!=1]; a,mg=conc(C)
    print(f"  {name}: conc={a:.4f}  |G(C)|={mg:.5f}  L>={mg*(1-a/7):.5f}")
print("\n=> if sup conc stays bounded < ~6.2 as Vmax grows, the band is a UNIFORM gap (quantization);")
print("   the AP (conc=7, regular circulant) is isolated above every covering family = the last bit.")
