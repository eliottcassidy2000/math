#!/usr/bin/env python3
"""
lrc14_family_completeness_search — mac-mini-2026-06-16-S2

The biggest gap in the inf L>0 program is FAMILY-COMPLETENESS: the numerics say
the interior-drop cores {1..13}\{j} ∪ {14m} are the extremizers, but a proof must
rule out ALL primitive 13-sets containing a multiple of 14. Here we test it
directly: exhaustively search 12-subsets of [1,N] (the "core") together with a
resonant stranger 14m, compute EXACT rational L, and check whether ANYTHING beats
the conjectured global min 1543/294294 (= {1..13}\{6} ∪ {98}).

If the min over this much larger family is still 1543/294294, that is strong
evidence the interior-drop cores are the true extremizers.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce

def danger_intervals(v):
    ivs=[]
    for k in range(0, v+1):
        lo=max(F(14*k-1,14*v),F(0)); hi=min(F(14*k+1,14*v),F(1))
        if lo<hi: ivs.append((lo,hi))
    return ivs

def union_measure(intervals):
    if not intervals: return F(0)
    iv=sorted(intervals); tot=F(0); cl,ch=iv[0]
    for lo,hi in iv[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: tot+=ch-cl; cl,ch=lo,hi
    return tot+(ch-cl)

_cache={}
def Lexact(S):
    key=tuple(sorted(set(S)))
    if key in _cache: return _cache[key]
    ivs=[]
    for v in key: ivs.extend(danger_intervals(v))
    r=F(1)-union_measure(ivs); _cache[key]=r; return r

GLOBAL_MIN=F(1543,294294)
strangers=[14*m for m in (1,2,3,4,5,7,11)]   # includes the resonant 98=14*7

def search(N):
    """All 12-subsets of [1,N]\\{multiples of 14}, each ∪ a resonant stranger."""
    pool=[v for v in range(1,N+1) if v%14!=0]
    best=[]   # (L, core, stranger)
    nbeat=0; ntot=0
    for core in combinations(pool,12):
        # primitivity of the whole set is automatic if 1 in core or gcd small; we
        # don't filter on primitivity (dilation only scales, doesn't lower L).
        for w in strangers:
            if w in core: continue
            S=list(core)+[w]
            L=Lexact(S); ntot+=1
            if L<GLOBAL_MIN: nbeat+=1
            best.append((L,core,w))
    best.sort(key=lambda t:t[0])
    return best[:8], nbeat, ntot

for N in (14,15,16):
    print("="*76)
    print(f"Search: 12-subsets of [1,{N}] (no mult of 14) ∪ resonant stranger 14m")
    print("="*76)
    top,nbeat,ntot=search(N)
    print(f"configs tested: {ntot};  configs with L < 1543/294294: {nbeat}")
    print(f"global min reference 1543/294294 = {float(GLOBAL_MIN):.9f}")
    print("smallest-L configurations found:")
    for L,core,w in top:
        drop=sorted(set(range(1,14))-set(core)) if set(core)<=set(range(1,14)) else None
        tag=f"  [= {{1..13}}\\{{{drop}}}]" if drop is not None else ""
        print(f"  L={float(L):.9f} = {L}  core={list(core)} +{w}{tag}")
    print()

print("="*76)
print("INTERPRETATION")
print("="*76)
print("If the global min over ALL these 12-subset families is STILL 1543/294294")
print("(attained at {1..13}\\{6} ∪ {98}), the interior-drop cores are confirmed")
print("extremal over this much larger search => strong evidence for completeness.")
print("Any config with L<1543/294294 would REFUTE the interior-drop conjecture.")
