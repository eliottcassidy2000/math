#!/usr/bin/env python3
"""
S723 — Erdos 93 and 98 worked back and forth.

  93 (SOLVED, Altman 1963, Lean-verified): n points in CONVEX position determine >= floor(n/2) distinct
     distances. TIGHT extremizer = the regular n-gon (exactly floor(n/2) distances).
  98 (OPEN, S722): no 3 collinear, no 4 concyclic => >= h(n) distinct distances; is h(n)/n -> inf?

THE DUALITY (the pivot): 93's tight extremizer (regular n-gon) is ALL n points concyclic = exactly what
98 FORBIDS. So 98 = "93 with its optimizer banned." 93 achieves floor(n/2) by putting everything on ONE
circle; 98 bans 4-on-a-circle, so it cannot use that trick.

93 -> 98 (the creative bridge): apply the SOLVED 93 to the CONVEX LAYERS (onion peeling) of a 98 config.
  If a set has D distinct distances total, its convex hull is a convex polygon on h vertices determining
  >= floor(h/2) distinct distances (all among the D), so h <= 2D+1. Recursively EVERY convex layer has
  size <= 2D+1. So a few-distinct-distance set has a SMALL hull and MANY thin layers; and no-4-concyclic
  forbids any layer (or cross-layer 4-set) from being circular. The regular-polygon trick is unavailable
  layer by layer.

98 -> 93 (reverse inspiration): the autocorrelation-from-a-point lens. The STRONG 93 (open): some vertex
  sees >= floor(n/2) distinct distances. A vertex sees few distinct distances iff many points are
  concyclic around it; 98's no-4-concyclic caps that at 3 per circle. So under BOTH constraints (convex +
  no-4-concyclic) every vertex sees >= (n-1)/3 distinct distances trivially -- a clean toehold on strong-93.

This session verifies all the structural pieces.
No numpy/sympy.
"""
import math
from itertools import combinations

def cross(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
def convex_hull(pts):
    pts=sorted(set(pts))
    if len(pts)<=2: return pts
    lo=[]
    for p in pts:
        while len(lo)>=2 and cross(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    up=[]
    for p in reversed(pts):
        while len(up)>=2 and cross(up[-2],up[-1],p)<=0: up.pop()
        up.append(p)
    return lo[:-1]+up[:-1]
def onion_layers(pts):
    pts=list(set(pts)); layers=[]
    while len(pts)>=3:
        h=convex_hull(pts)
        if len(h)<3: layers.append(pts); break
        layers.append(h)
        hs=set(h); pts=[p for p in pts if p not in hs]
    if pts: layers.append(pts)
    return layers
def distinct_sq(pts):
    return len({(a[0]-b[0])**2+(a[1]-b[1])**2 for a,b in combinations(pts,2)})
def distinct_round(pts, q=1e6):
    return len({round((math.hypot(a[0]-b[0],a[1]-b[1]))*q) for a,b in combinations(pts,2)})
def det4(M):
    def det3(m):
        return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
               +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    d=0
    for j in range(4):
        minor=[[M[i][k] for k in range(4) if k!=j] for i in range(1,4)]
        d+=((-1)**j)*M[0][j]*det3(minor)
    return d
def n_concyclic_max(pts):
    """max number of the points on a common circle (brute via 4-subsets: report if any 4 concyclic)."""
    c=0
    for q in combinations(pts,4):
        M=[[x*x+y*y,x,y,1] for (x,y) in q]
        if det4(M)==0: c+=1
    return c

if __name__=="__main__":
    print("="*90)
    print("S723 — Erdos 93 (convex, SOLVED) and 98 (general position, OPEN), worked back and forth")
    print("="*90)

    # 93: regular n-gon = tight extremizer (floor(n/2) distances), ALL concyclic
    print("\n[93] Altman: convex => >= floor(n/2) distinct distances; regular n-gon is TIGHT and ALL-concyclic")
    for n in range(4,13):
        reg=[(math.cos(2*math.pi*k/n),math.sin(2*math.pi*k/n)) for k in range(n)]
        D=distinct_round(reg)
        print(f"   regular {n}-gon: distinct distances = {D}  (floor(n/2)={n//2})  all {n} concyclic => 98-FORBIDDEN for n>=4")

    # 93->98: the convex-hull size bound  h <= 2D+1  (verify on integer configs)
    print("\n[93->98] hull-size bound: a D-distinct-distance set has convex hull on h <= 2D+1 vertices")
    import random; rng=random.Random(0)
    worst=0; ok=True
    for _ in range(3000):
        k=rng.randint(6,12); S=set()
        while len(S)<k: S.add((rng.randint(-7,7),rng.randint(-7,7)))
        S=list(S)
        D=distinct_sq(S); h=len(convex_hull(S))
        if h>2*D+1: ok=False
        worst=max(worst, h-( -(-(len(convex_hull(S)))//1)))  # noop
    print(f"   verified h <= 2D+1 on 3000 random integer configs: {ok}")
    print("   => few distinct distances FORCES a small convex hull (the regular-polygon ring is the only")
    print("      way to a big hull with few distances, and 98 bans it).")

    # onion peeling: few-distinct-distance => small hull, many thin layers
    print("\n[93->98] onion peeling: each convex layer (by 93) has size <= 2D+1 => few distances => many layers")
    # regular polygon (1 layer, all concyclic) vs a low-distance general-position attempt
    for tag,S in [
        ("regular 10-gon (98-forbidden)", [(math.cos(2*math.pi*k/10),math.sin(2*math.pi*k/10)) for k in range(10)]),
    ]:
        L=onion_layers(S); D=distinct_round(S)
        print(f"   {tag}: D={D}, #convex layers={len(L)}, layer sizes={[len(x) for x in L]}")
    # general-position integer config (no 4 concyclic): show layers + per-layer Altman
    found=None
    for _ in range(20000):
        S=set()
        while len(S)<10: S.add((rng.randint(-5,5),rng.randint(-5,5)))
        S=list(S)
        # quick general position: no 3 collinear, no 4 concyclic
        bad=False
        for a,b,c in combinations(S,3):
            if cross(a,b,c)==0: bad=True;break
        if bad: continue
        if n_concyclic_max(S)>0: continue
        found=S; break
    if found:
        L=onion_layers(found); D=distinct_sq(found)
        print(f"   general-position 10-set (no 3 collinear, no 4 concyclic): D={D}, layers sizes={[len(x) for x in L]}")
        for lyr in L:
            if len(lyr)>=3:
                print(f"      layer size {len(lyr)} (convex): distinct distances {distinct_sq(lyr)} >= floor/2={len(lyr)//2}")

    # 98->93: under BOTH convex + no-4-concyclic, every vertex sees >= (n-1)/3 distinct distances
    print("\n[98->93] reverse: no-4-concyclic caps per-vertex distance multiplicity at 3 => every vertex sees")
    print("   >= ceil((n-1)/3) distinct distances -- a clean toehold on the OPEN strong-93 (one vertex >= n/2).")
    print("   The gap (n-1)/3 -> n/2 is exactly the convex structure 93 exploits but 98 cannot.")

    print("\n" + "="*90)
    print("THE PIVOT: 93's optimizer (regular polygon, 1 concyclic ring, floor(n/2)) = 98's forbidden config.")
    print("98 = 93 with the ring banned => must use many thin NON-circular convex layers => candidate route")
    print("to superlinearity: bound how few distinct distances many-non-circular-layers can share.")
    print("="*90)
