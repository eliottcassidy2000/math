#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Map the LRC(14) TIGHT LOCUS (primitive 13-sets with L=0) and the TRUE inf of the lonely
measure over ALL primitive 13-sets — testing whether the repo's inf≈0.0052 (restricted to
multiple-of-14 strangers) is the truth, or whether perturbations of SPORADIC tight configs
go lower (already found {1..11,13,36}=1/1260≈0.000794) and possibly →0 (which would kill C'(14)).
kind-pasteur-2026-06-16-S6.  L(S)=meas{τ:||vτ||>1/14 ∀v} computed EXACTLY (rational arc sweep).
"""
import sys; sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools

def danger_arcs(v):
    w=Fr(1,14*v); arcs=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: arcs.append((Fr(0),hi)); arcs.append((1+lo,Fr(1)))
        elif hi>1: arcs.append((lo,Fr(1))); arcs.append((Fr(0),hi-1))
        else: arcs.append((lo,hi))
    return arcs
def union_measure(arcs):
    arcs=sorted((a,b) for a,b in arcs if b>a); tot=Fr(0); cl=ch=None
    for a,b in arcs:
        if ch is None: cl,ch=a,b
        elif a<=ch: ch=max(ch,b)
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return tot
def L_exact(S):
    arcs=[]
    for v in S: arcs+=danger_arcs(v)
    return 1-union_measure(arcs)
def primitive(S): return reduce(gcd,S)==1

if __name__=="__main__":
    # ---- (A) the TIGHT LOCUS via 1- and 2-element replacement of {1..13} ----
    print("=== (A) tight locus (L=0) among 1- and 2-element replacements of {1..13}, w<=120 ===")
    tights=[]
    for j in range(1,14):
        core=[x for x in range(1,14) if x!=j]
        for w in range(14,121):
            if w in core: continue
            S=tuple(sorted(core+[w]))
            if not primitive(S): continue
            if L_exact(S)==0: tights.append(S)
    print(f"   1-drop tight configs (L=0): {len(tights)} found")
    for S in tights[:20]: print(f"     {S}")
    # 2-drop tights (smaller w range)
    t2=0
    for j1,j2 in itertools.combinations(range(1,14),2):
        core=[x for x in range(1,14) if x not in (j1,j2)]
        for w1 in range(14,60):
            for w2 in range(w1+1,80):
                S=tuple(sorted(core+[w1,w2]))
                if len(set(S))!=13 or not primitive(S): continue
                if L_exact(S)==0: t2+=1
    print(f"   2-drop tight configs (L=0, w<=79): {t2} found")

    # ---- (B) TRUE inf: min positive L over 1-drop (large w) and 2-drop ----
    print("\n=== (B) hunting the true min positive L (is it bounded, or →0?) ===")
    best=(Fr(1),None)
    for j in range(1,14):
        core=[x for x in range(1,14) if x!=j]
        for w in range(14,400):
            if w in core: continue
            S=tuple(sorted(core+[w]))
            if not primitive(S): continue
            Lv=L_exact(S)
            if 0<Lv<best[0]: best=(Lv,S)
    print(f"   1-drop, w<400: min positive L = {float(best[0]):.8f} = {best[0]} at {best[1]}")
    b1=best
    for j1,j2 in itertools.combinations(range(1,14),2):
        core=[x for x in range(1,14) if x not in (j1,j2)]
        for w1 in range(14,80):
            for w2 in range(w1+1,120):
                S=tuple(sorted(core+[w1,w2]))
                if len(set(S))!=13 or not primitive(S): continue
                Lv=L_exact(S)
                if 0<Lv<best[0]: best=(Lv,S)
    print(f"   incl 2-drop (w<=119): min positive L = {float(best[0]):.8f} = {best[0]} at {best[1]}")

    # ---- (C) perturb each found tight config minimally; smallest opened L ----
    print("\n=== (C) minimal perturbations of tight configs → smallest opened L ===")
    pbest=(Fr(1),None)
    for T in tights:
        for idx in range(13):
            for delta in (-3,-2,-1,1,2,3,6,7,12,14):
                S=list(T); S[idx]=S[idx]+delta
                if len(set(S))!=13 or min(S)<1 or not primitive(tuple(S)): continue
                Lv=L_exact(tuple(sorted(S)))
                if 0<Lv<pbest[0]: pbest=(Lv,tuple(sorted(S)))
    print(f"   smallest L from minimal tight-perturbations: {float(pbest[0]):.8f} = {pbest[0]} at {pbest[1]}")
    gmin=min(best,b1,pbest,key=lambda x:x[0])
    print(f"\n=== GLOBAL min positive L found this session: {float(gmin[0]):.8f} = {gmin[0]} at {gmin[1]} ===")
    print(f"   (repo's stated inf was 0.0052; 1/1260={float(Fr(1,1260)):.8f}; is the true min below 1/1260?)")
