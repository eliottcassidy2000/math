#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CORE -> CIRCULANT-GRAPH MAP + DESCENT GAP-PROFILES (klein-S21, Part C+D).

Each core O has autocorrelation a(d)=#{x:x,x+d in O}; the gap = lambda_min of the circulant [a(i-j)]
= lambda_min of a CIRCULANT GRAPH (signless-Laplacian-type). Map each gap class to its graph type, and
order cores by 'spread'. Then the DESCENT GAP-PROFILE of a covering = the multiset of gaps along its
descent chain (a new fine trackable).
"""
import math, cmath, itertools
P=7; W=cmath.exp(2j*math.pi/P)
def autocorr(O):
    Os=set(O); return [sum(1 for x in Os if (x+d)%P in Os) for d in range(P)]
def gap(O):
    sp=[abs(sum(W**((k*x)%P) for x in set(O)))**2 for k in range(P)]
    return min(sp[k] for k in range(1,P))

print("="*86); print(" CORE -> CIRCULANT-GRAPH MAP: each gap class is a named circulant graph"); print("="*86)
examples=[("full Z_7","7J (rank-1)",set(range(7))),
          ("doublet {0,1}","2I+A(C_7) = Q(C_7) signless Laplacian of the 7-CYCLE",{0,1}),
          ("arc {0,1,2}","interval autocorr [3,2,1,0,0,1,2]",{0,1,2}),
          ("singleton {0}","I (empty graph)",{0}),
          ("diff-set/QR {1,2,4}","2I+(J-I)=I+J: COMPLETE graph K_7 (Paley)",{1,2,4})]
for nm,desc,O in examples:
    print(f"   {nm:<16} a(d)={autocorr(O)}  gap={gap(O):.4f}  <- {desc}")
print("\n   => LOW gap = CONCENTRATED core (doublet 0.198 < arc 0.308); HIGH gap = SPREAD (diff-set=2).")
print("      the gap ORDERS cores by spread; the binding obstruction is the most-concentrated (doublet).")

print("\n"+"="*86); print(" SPREAD ORDERING (gap as a concentration index), all Z_7 core shapes:"); print("="*86)
shapes=[("doublet (2 adj)",{0,1}),("pair-skip1 (2, d=2)",{0,2}),("arc-3 {0,1,2}",{0,1,2}),
        ("V-3 {0,1,3}",{0,1,3}),("diff-set {1,2,4}",{1,2,4}),("singleton",{0}),
        ("arc-4 {0,1,2,3}",{0,1,2,3}),("co-doublet (5)",{0,1,2,3,4})]
for nm,O in sorted(shapes,key=lambda t:gap(t[1])):
    print(f"   gap={gap(O):.4f}  {nm:<22} |O|={len(O)}")

print("\n"+"="*86); print(" DESCENT GAP-PROFILE (new trackable): the multiset of apex gaps along the descent chain"); print("="*86)
def descend(S):
    cores=[]; cur=sorted(set(S)); g=0
    while cur and g<30:
        O=[v for v in cur if v%2==1]; E=[v for v in cur if v%2==0]; cores.append(O)
        if not E: break
        cur=sorted(set(v//2 for v in E)); g+=1
    return cores
covs={"consec{1..13}":list(range(1,14)),"tightest{1..12,182}":list(range(1,13))+[182],
      "skip12{1..11,13,84}":list(range(1,12))+[13,84],"binding{1..13}\\7":[x for x in range(1,14) if x!=7],
      "consec{1..14}":list(range(1,15)),"consec{1..7}":list(range(1,8))}
for nm,S in covs.items():
    prof=[]
    for O in descend(S):
        if O: prof.append(round(gap(O),3))
    minpos=min([g for g in prof if g>1e-9], default=None)
    print(f"   {nm:<22} gap-profile {prof}   min-positive={minpos}  (depth {len(prof)})")
print("   => every chain's min-positive gap is the DOUBLET 0.198 (or a 0.308) -- the binding per-level atom;")
print("      a gap of 0.0 in the profile = the chain passes through the full-Z_7 cusp (mod-7 covering at that level).")
