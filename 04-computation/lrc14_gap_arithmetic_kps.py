#!/usr/bin/env python3
"""
lrc14_gap_arithmetic_kps  (kind-pasteur, PROVE side)

DERIVE the arithmetic of the 1/1260 champion gap from first principles.

The AP {1..13} drop-12 exclusive gap G_12 consists of 4 arcs. We show each arc
is centered at a 12-resonance τ = (14k±1)/(14·12) endpoint configuration, and
explain exactly why w=36 leaves residual length 1/2520 in the two central arcs.

KEY ARITHMETIC FACT we want to expose:
  - A residual lonely arc has length 1/2520 = 1/(14 · 180) = 1/(14·lcm-related).
  - 1260 = 14 · 90 = 2^2·3^2·5·7. 2520 = 2·1260.
  - The residual arc lives where 12's danger arc (half-width 1/(14·12)=1/168) was the
    ONLY cover, and 36 (half-width 1/(14·36)=1/504) only re-covers part because
    36's danger centers k/36 don't all align with 12's centers k/12.

We make this explicit: list 12's danger-arc centers in G_12, and for each, the
nearest 36-center and the residual.
"""
from fractions import Fraction as F

def danger_arcs_raw(v):
    """list of (center=k/v, lo, hi) including wraps, k=0..v"""
    w=F(1,14*v); out=[]
    for k in range(v+1):
        c=F(k,v); out.append((c,c-w,c+w))
    return out
def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(a,b) for a,b in out if b>a]
def union(arcs):
    arcs=sorted((a,b) for a,b in arcs if b>a)
    if not arcs: return []
    res=[]; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=max(ch,hi)
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch)); return res
def total(a): return sum(b-x for x,b in union(a))
def complement(arcs):
    u=union(arcs); res=[]; prev=F(0)
    for lo,hi in u:
        if lo>prev: res.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: res.append((prev,F(1)))
    return res
def gap_of_drop(C,e):
    cov=[]
    for v in C:
        if v!=e: cov+=danger(v)
    return complement(cov)

C=list(range(1,14))
G12=gap_of_drop(C,12)
print("G_12 arcs (the points ONLY speed 12 protected in AP):")
for lo,hi in G12:
    # which 12-center is this near?
    c12=min((abs(lo-F(k,12)),F(k,12)) for k in range(13))[1]
    print(f"  [{lo},{hi}] len={hi-lo}  near 12-center k/12={c12}={float(c12):.4f}")

print("\nWhy these 4 arcs? 12's danger centers are k/12. The gap is the part of")
print("12's danger arc [k/12-1/168, k/12+1/168] NOT covered by speeds 1..11,13.")
print("12's arc half-width = 1/168. Listing the 4 gap arcs vs the 12-arc they sit in:\n")

# 36's danger arcs
D36=danger(36)
print("Now ADD w=36 (half-width 1/(14*36)=1/504). Residual = G_12 minus D_36:")
resid=complement([a for a in danger(36)] + [a for v in C if v!=12 for a in danger(v)])
# residual lonely = complement of all-arcs of perturbed config
S=[x for x in C if x!=12]+[36]
allarcs=[]
for v in S: allarcs+=danger(v)
lonely=complement(allarcs)
print(f"  total residual lonely L = {total(allarcs) and (F(1)-total(allarcs))}")
print(f"  lonely arcs:")
for lo,hi in lonely:
    # express endpoints in terms of /504 and /2520
    print(f"    [{lo},{hi}]  len={hi-lo}")
    # decode endpoints
    for label,x in [("lo",lo),("hi",hi)]:
        # find which speed's danger boundary this is
        for v in S:
            wv=F(1,14*v)
            for k in range(v+1):
                if x==F(k,v)-wv or x==F(k,v)+wv:
                    sgn = '-' if x==F(k,v)-wv else '+'
                    print(f"        {label}={x} = {k}/{v} {sgn} 1/{14*v}  (speed {v}, center {k}/{v})")

print("\nINTERPRETATION:")
print(" The two residual arcs are bounded on one side by speed-13's danger arc")
print(" and on the other by speed-36's danger arc. The gap between consecutive")
print(" danger centers of 13 and 36 near tau=29/70 is exactly 1/2520 wide.")
print()
# explicit: the left residual [29/70, 209/504]
lo,hi=lonely[0]
print(f" Left residual [{lo},{hi}]: lo=29/70, hi=209/504.")
print(f"   29/70 = boundary of speed-? ; 209/504=boundary of speed-36 (209/504=?)")
# 29/70 as k/v - 1/14v
for v in S:
    wv=F(1,14*v)
    for k in range(v+1):
        if F(29,70)==F(k,v)-wv: print(f"   29/70 = {k}/{v} - 1/{14*v}  (right end of speed-{v} arc just left)")
        if F(29,70)==F(k,v)+wv: print(f"   29/70 = {k}/{v} + 1/{14*v}  (right end of speed-{v} arc)")
for v in S:
    wv=F(1,14*v)
    for k in range(v+1):
        if F(209,504)==F(k,v)-wv: print(f"   209/504 = {k}/{v} - 1/{14*v}  (left end of speed-{v} arc just right)")

print()
print("LCM analysis of champion config {1..11,13,36}:")
import math
S2=[x for x in C if x!=12]+[36]
l=1
for v in S2: l=l*v//math.gcd(l,v)
print(f"  lcm = {l}, 14*lcm = {14*l}, 1/(14*lcm)={F(1,14*l)}={float(F(1,14*l)):.3e}")
print(f"  L=1/1260; 1260 divides 14*lcm? {14*l} / 1260 = {F(14*l,1260)}")
print(f"  quantization: L in (1/(14*lcm))Z; 1/1260 = {F(14*l,1260)} units of 1/(14*lcm)")
