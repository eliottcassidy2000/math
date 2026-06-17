#!/usr/bin/env python3
"""
LRC(14) DISPROVE - near-miss combination search.
The obstruction analysis found single perturbations giving M=2/27=0.07407 (very
close to 1/14=0.07143): {1..13} with 13->26 or 10->20. Also the doubling 12->24
stays EXACTLY at 1/14. We now COMBINE near-miss moves and search their joint
neighborhood to see if any combination dips below 1/14.

SEED FROM PROVE: the prove route shows tight minimizers are pinned at tau=k/14 by
extreme-speed pairs ({1,13},{5,9},{3,11} bind 1/14,3/14,5/14). A counterexample
must unpin ALL THREE plateaus at once. We target moves that touch the binding
speeds {1,3,5,9,11,13} simultaneously.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations, product
import random

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def M(S):
    b=F(0); at=None
    for t in cand(S):
        v=g(S,t)
        if v>b: b=v; at=t
    return b,at
def is_prim(S): return reduce(gcd,S)==1

THRESH=F(1,14)
ce=[]; best=[None]; seen=set()
def consider(S):
    S=tuple(sorted(set(int(x) for x in S)))
    if len(S)!=13 or any(v<=0 for v in S): return
    if not is_prim(S): return
    if S in seen: return
    seen.add(S)
    m,at=M(S)
    if best[0] is None or m<best[0][0]: best[0]=(m,S,at)
    if m<THRESH: ce.append((m,S,at))

base=list(range(1,14))
binders=[1,3,5,9,11,13]  # the speeds pinning 1/14,3/14,5/14

# (1) target the binding speeds: replace each binder with a near-multiple,
# try all pairs/triples of binder replacements.
print("(1) replace binding speeds (pins of 1/14,3/14,5/14) jointly", flush=True)
# candidate replacements for a binder v: 2v, 2v+-1, v+14, v+-7
def repls(v):
    out=set()
    for r in [2*v,2*v-1,2*v+1,v+14,v+7,v+13,v+15,3*v]:
        if r>0: out.add(r)
    return out
# single binder replacement (full scan)
for b in binders:
    idx=base.index(b)
    for r in range(14,40):
        consider(base[:idx]+[r]+base[idx+1:])
# pairs of binder replacements
for (b1,b2) in combinations(binders,2):
    i1,i2=base.index(b1),base.index(b2)
    for r1 in repls(b1):
        for r2 in repls(b2):
            S=base[:]; S[i1]=r1; S[i2]=r2
            consider(S)
# triples
for (b1,b2,b3) in combinations(binders,3):
    i1,i2,i3=base.index(b1),base.index(b2),base.index(b3)
    for r1 in list(repls(b1))[:4]:
        for r2 in list(repls(b2))[:4]:
            for r3 in list(repls(b3))[:4]:
                S=base[:]; S[i1]=r1; S[i2]=r2; S[i3]=r3
                consider(S)
print(f"  ce={len(ce)} best={best[0][0]}={float(best[0][0]):.6f}", flush=True)

# (2) climb from near-miss champions 13->26 and 10->20: local search
print("(2) local search from near-miss champions (M=2/27)", flush=True)
champs=[tuple(sorted(set(base[:12]+[26]))),    # 13->26
        tuple(sorted(set(base[:9]+base[10:]+[20])))]  # 10->20
random.seed(2026)
frontier=set()
for c in champs:
    if len(c)==13: frontier.add(c)
for _ in range(50):
    newf=set()
    for S in list(frontier)[:200]:
        L=list(S)
        for i in range(13):
            for delta in (-2,-1,1,2,7,-7,13,14):
                v=L[i]+delta
                if v<=0 or v in L: continue
                cand_set=tuple(sorted(set(L[:i]+[v]+L[i+1:])))
                if len(cand_set)==13 and is_prim(cand_set):
                    m,at=M(cand_set)
                    if best[0] is None or m<best[0][0]: best[0]=(m,cand_set,at)
                    if m<THRESH: ce.append((m,cand_set,at))
                    # keep promising (near 1/14) on frontier
                    if m < F(2,27):
                        newf.add(cand_set)
    if not newf: break
    frontier=newf
print(f"  ce={len(ce)} best={best[0][0]}={float(best[0][0]):.6f}", flush=True)

# (3) scaled tight AP with a few unit-shifts: c*{1..13} +- small, primitive
print("(3) scaled-and-shifted tight cores", flush=True)
for c in range(1,8):
    core=[c*k for k in range(1,14)]
    # shift a few elements by +-1,+-2 to keep primitive
    for shifts in product([-2,-1,0,1,2],repeat=3):
        for trip in combinations(range(13),3):
            S=core[:]
            ok=True
            for idx,sh in zip(trip,shifts):
                S[idx]=core[idx]+sh
                if S[idx]<=0: ok=False
            if ok: consider(S)
        break  # only first triple-position set per c to bound cost; rely on (1)/(2)
print(f"  ce={len(ce)} best={best[0][0]}={float(best[0][0]):.6f}", flush=True)

print("="*70, flush=True)
print(f"Counterexamples (M<1/14): {len(ce)}", flush=True)
ce.sort()
for m,S,at in ce[:25]:
    print(f"  M={m}={float(m):.6f} S={S} at {at}", flush=True)
m,S,at=best[0]
print(f"\nSMALLEST M: {m}={float(m):.8f}  S={S} tau={at}", flush=True)
print(f"  M-1/14={m-THRESH}={float(m-THRESH):.8f}; M>=1/14? {m>=THRESH}", flush=True)
