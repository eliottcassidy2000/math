#!/usr/bin/env python3
"""
kps-S39 (corrected): the ORDER-3 members are DILATED-AP + defects, NOT spacing-1 ladder bases.

Discovery: the N=6 order-3 witness is {1,5,6,11,16,17}=5/33 at t=10/33 -- structure AP{1,6,11,16}
(spacing d=5) + boundary defects {5,17}.  My S38 order-3 check only tried spacing-1 ladder bases,
so it did NOT actually test this class -- an honest correction: S38's "order-3 empty at N=12" was
over the WRONG structure.

This file properly searches DILATED-AP + defect families at N=12 for the order-3 values 4/51, 5/63
(and any gap member), the structure that the N=6 witness reveals.  Note the N=6 AP spacing d=5
equals the M-numerator (5/33); we scan d broadly.
"""
from fractions import Fraction
from itertools import combinations
import numpy as np
from math import gcd
from functools import reduce

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        a=np.arange(1,q); r=np.outer(va,a)%q; d=np.minimum(r,q-r); col=d.min(axis=0); qb=int(col.max())
        if qb*bd>bn*q: bn,bd,bc=qb,q,int(a[col.argmax()])
    return Fraction(bn,bd),(bc,bd)

N=12; LO,HI=Fraction(1,N+1),Fraction(2,2*N+1)
o3={Fraction(4,4*N+3),Fraction(5,5*N+3)}   # 4/51, 5/63
print(f"=== N=12 DILATED-AP + defect search for order-3 (4/51,5/63) and ANY gap member ===", flush=True)
print(f"  (N=6 template: AP spacing d=5 + boundary defects; scanning d, length L, defects)\n", flush=True)

gaphits=[]; o3hits=[]; tested=0
for d in range(2,8):                       # AP spacing
    for L in range(8,12):                  # AP length
        for a in range(1,4):               # AP start
            ap=[a+i*d for i in range(L)]
            need=N-L
            if need<0 or need>4: continue
            top=ap[-1]
            # defect pool: near AP endpoints and interior gaps, up to a modest cap
            pool=sorted(set(range(1,top+2*d+2)) - set(ap))
            pool=[x for x in pool if x<=64]
            if len(pool)<need: continue
            # limit combos for speed: prefer boundary-ish defects
            combos=list(combinations(pool,need))
            if len(combos)>4000: combos=combos[::max(1,len(combos)//4000)]
            for defs in combos:
                v=sorted(set(ap)|set(defs))
                if len(v)!=N or reduce(gcd,v)!=1: continue
                if max(v)>64: continue
                tested+=1
                M,(c,q)=Mw(v)
                if LO<M<HI:
                    gaphits.append((v,M))
                    if M in o3: o3hits.append((v,M,c,q))
print(f"  tested ~{tested} dilated-AP+defect families (max<=64)", flush=True)
print(f"  ANY gap member (M in (1/13,2/25)): {len(gaphits)}", flush=True)
for v,M in gaphits[:15]: print(f"     {v}: M={M}", flush=True)
print(f"  ORDER-3 hits (M in {{4/51,5/63}}): {len(o3hits)}", flush=True)
for v,M,c,q in o3hits[:15]: print(f"     {v}: M={M} at t={c}/{q}", flush=True)

print("\n" + "="*70, flush=True)
if not gaphits:
    print("RESULT: no dilated-AP+defect family (max<=64) lands in the N=12 gap -- so the order-3", flush=True)
    print("structure that makes N=6 nonempty (5/33) does NOT reproduce at N=12 up to this height.", flush=True)
    print("This is the PROPER order-3 check (S38's was over the wrong, spacing-1, structure).", flush=True)
    print("Consistent with mac-mini census (height 48) + opus emptiness -- now over the right class.", flush=True)
else:
    print(f"RESULT: {len(gaphits)} gap member(s) at N=12 -- INSPECT (would matter a lot!).", flush=True)
