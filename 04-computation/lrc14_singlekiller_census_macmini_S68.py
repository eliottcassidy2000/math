#!/usr/bin/env python3
"""mac-mini-S68 (lean): confirm no PRIMITIVE COVERING single-killer config beats 14/183.
Early-exit: a config CLEARS if some q<=Qmax gives min-residue-dist >= 14q/183. Focus on the
near-deep-well cores (one element of {1..12} pushed out) + dilated cores -- the candidates that
could approach the floor."""
from fractions import Fraction as F
from math import gcd
target=F(14,183)
def clears_or_M(S,Qmax):
    """return (clears>=target?, best_M_seen). Early-exit on first q with best_q>=target."""
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=q; ok=True
            for v in S:
                r=(a*v)%q; d=r if r<=q-r else q-r
                if d<mind: mind=d
                if F(d,q)<target: ok=False; break
            if ok and mind>0: return True,F(mind,q)
            if mind>0 and F(mind,q)>best: best=F(mind,q)
    return False,best
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

cores=[]
for drop in range(1,13):
    for add in range(13,60):
        cores.append(sorted(set([x for x in range(1,13) if x!=drop]+[add])))
for c in range(1,8): cores.append([c*k for k in range(1,13)])
# two-swap near-deep-well
for d1 in range(1,13):
    for d2 in range(d1+1,13):
        for add in [ (13,14),(14,15),(13,26),(24,26) ]:
            cores.append(sorted(set([x for x in range(1,13) if x not in (d1,d2)])|set(add)))
killers=[182,364,546,910,1274]
tested=0; below=[]; nonclear=[]; mins=[]
seen=set()
for core in cores:
    if len(core)!=12: continue
    tc=tuple(core)
    if tc in seen: continue
    seen.add(tc)
    for killer in killers:
        S=sorted(set(core+[killer]))
        if len(S)!=13 or not prim(S) or not is_cov(S): continue
        tested+=1
        cl,best=clears_or_M(S,220)
        if not cl:
            cl2,best2=clears_or_M(S,2*max(S)+5)
            if not cl2: below.append((best2,S))
            else: mins.append((best2,S))
        else: mins.append((best,S))
mins.sort()
print(f"tested {tested} primitive covering single-killer configs (near-deep-well + dilated).")
print(f"floor 14/183={float(target):.6f}")
print("below floor:", "NONE" if not below else below[:6])
if mins:
    print(f"global min M = {float(mins[0][0]):.6f}={mins[0][0]} at S={mins[0][1]}")
    print("smallest few:")
    for M,S in mins[:8]:
        tag=" <== DEEP WELL" if S==sorted([*range(1,13),182]) else ""
        print(f"   M={float(M):.6f}={M} S={S}{tag}")
