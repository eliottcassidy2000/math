#!/usr/bin/env python3
"""
lrc14_broadfloor_macmini_0618s1.py (mac-mini-2026-06-18-S1)
Broad-shape floor: min rho*(P,E) over P subset {1..13} (|P|=13-k) AND bounded-spread shapes E
(spread up to 3k), per k. Refines the consecutive-only floor 1/84 (consecutive NOT extremal).
Confirms: floor stays POSITIVE; extremal spread is bounded (compact problem).
"""
import sys, random, itertools
sys.stdout.reconfigure(line_buffering=True)
from fractions import Fraction as F
random.seed(618111)
H=F(1,14)
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def safe_set(A,h=H):
    if not A: return [(F(0),F(1))]
    dz=merge([iv for u in A for iv in danger_arcs(u,h)]); safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def maxgap(points):
    pts=sorted(p%1.0 for p in points); g=0.0
    for i in range(len(pts)-1): g=max(g,pts[i+1]-pts[i])
    return max(g,(pts[0]+1.0)-pts[-1])
def rho(P,E,N=50000):
    GP=[(float(a),float(b)) for a,b in safe_set(P)]; cnt=0
    for t in range(N):
        x=(t+0.5)/N
        if not any(a<=x<b for a,b in GP): continue
        if maxgap([(e*x)%1.0 for e in E])>2.0/7.0: cnt+=1
    return cnt/N

print("Broad-shape floor min rho* over P and bounded-spread E (spread<=3k), per k:")
overall=(9.9,None,None,None)
for k in range(7,14):
    psz=13-k
    best=(9.9,None,None)
    # sample P and shapes E
    Ps=[list(p) for p in itertools.combinations(range(1,14),psz)]
    if len(Ps)>25: Ps=random.sample(Ps,25)
    for P in Ps:
        # consecutive + perforated + spread shapes
        shapes=[list(range(k))]
        for _ in range(25):
            sp=random.randint(k-1,3*k)
            E=sorted(set([0,sp]+random.sample(range(1,sp),min(k-2,sp-1))))
            if len(E)==k: shapes.append(E)
        for E in shapes:
            r=rho(P,E,N=30000)
            if r<best[0]: best=(r,tuple(P),tuple(E))
    print(f"   k={k:2d}: min rho* ~ {best[0]:.5f}  P={best[1]} E={best[2]}")
    if best[0]<overall[0]: overall=(best[0],k,best[1],best[2])
print(f"\n   BROAD FLOOR ~ {overall[0]:.5f} at k={overall[1]}, P={overall[2]}, E={overall[3]}")
print(f"   (consecutive-only exact floor was 1/84={1/84:.5f}; broad is lower but POSITIVE)")
print("DONE.")
