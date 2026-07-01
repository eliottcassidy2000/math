#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Does (I(Omega,x), d=det(I+S)) stay INJECTIVE on n=7 tournament classes? (or does a cospectral family reappear
and force a THIRD reconstruction coordinate?)

kind-pasteur-2026-07-01-S15. At n=6 the fingerprint (I(Omega,x), d) separated all classes (the OCF-cospectral
twin {4,6} was d-visible). Test n=7: sample tournament classes, compute (I(Omega,x), d [, skew spectrum cpS]),
and report any (I(Omega,x), d)-collision between NON-isomorphic classes.
"""
import sys, itertools, random
import numpy as np
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
n=7; P=list(itertools.permutations(range(n)))

def canon(A):
    best=None
    for p in P:
        s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or s<best: best=s
    return best
def rand_tour(rng):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<0.5: A[i][j]=1
            else: A[j][i]=1
    return A
def odd_cycle_verts(A):
    out=[]
    for L in range(3,n+1,2):
        for sub in itertools.combinations(range(n),L):
            s0=sub[0]
            for pr in itertools.permutations(sub[1:]):
                seq=(s0,)+pr
                if all(A[seq[t]][seq[(t+1)%L]] for t in range(L)):
                    out.append((min(tuple(seq[k:]+seq[:k]) for k in range(L)), frozenset(sub)))
    seen=set(); vs=[]
    for key,fs in out:
        if key not in seen: seen.add(key); vs.append(fs)
    return vs
def fingerprint(A):
    vs=odd_cycle_verts(A); i1=len(vs)
    i2=sum(1 for a,b in itertools.combinations(range(i1),2) if vs[a].isdisjoint(vs[b]))
    i3=0
    for a,b,c in itertools.combinations(range(i1),3):
        if vs[a].isdisjoint(vs[b]) and vs[a].isdisjoint(vs[c]) and vs[b].isdisjoint(vs[c]): i3+=1
    M=np.array(A); S=M-M.T
    d=int(round(np.linalg.det(np.eye(n)+S)))
    cpS=tuple(int(round(c)) for c in np.poly(S))
    ipoly=(i1,i2,i3)
    return ipoly,d,cpS

rng=random.Random(7)
NS=int(sys.argv[1]) if len(sys.argv)>1 else 1400
byCanon={}   # canon -> (ipoly,d,cpS)
for it in range(NS):
    A=rand_tour(rng); c=canon(A)
    if c not in byCanon: byCanon[c]=fingerprint(A)
classes=list(byCanon.items())
print(f"n=7: sampled {NS} random tournaments -> {len(classes)} distinct iso classes found (of 456, "
      f"{100*len(classes)/456:.0f}% coverage)")
# check (ipoly,d) collisions among distinct classes
from collections import defaultdict
byID=defaultdict(list); byIDd=defaultdict(list)
for c,(ip,d,cpS) in classes:
    byID[ip].append(c)           # OCF poly only
    byIDd[(ip,d)].append(c)      # OCF poly + determinant
ocf_coll=[(k,v) for k,v in byID.items() if len(v)>1]
ocfd_coll=[(k,v) for k,v in byIDd.items() if len(v)>1]
print(f"  (I(Omega,x)) alone: {len(ocf_coll)} cospectral groups (>=2 classes share the OCF poly) "
      f"covering {sum(len(v) for k,v in ocf_coll)} classes")
print(f"  (I(Omega,x), d): {len(ocfd_coll)} collision groups "
      f"=> (I(Omega,x),d) INJECTIVE on the sample? {len(ocfd_coll)==0}")
if ocfd_coll:
    print("   FIRST (I(Omega,x),d)-cospectral collisions (need a 3rd coordinate):")
    for (ip,d),v in ocfd_coll[:5]:
        print(f"     ipoly(i1,i2,i3)={ip}, d={d}: {len(v)} classes; skew-spectra distinct? "
              f"{len(set(byCanon[c][2] for c in v))==len(v)}")
    # does adding cpS (skew spectrum) resolve them?
    byIDdS=defaultdict(list)
    for c,(ip,d,cpS) in classes: byIDdS[(ip,d,cpS)].append(c)
    dS_coll=[v for k,v in byIDdS.items() if len(v)>1]
    print(f"   adding skew-spectrum cpS: {len(dS_coll)} residual collision groups "
          f"=> (I(Omega,x),d,cpS) injective on sample? {len(dS_coll)==0}")
else:
    print("   => on this sample, the (OCF, determinant) pair is a COMPLETE fingerprint at n=7 (no 3rd coord needed yet).")
print("DONE.")
