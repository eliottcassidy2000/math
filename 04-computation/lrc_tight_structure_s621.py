#!/usr/bin/env python3
"""
S621 part 2 — STRUCTURE of the exactly-tight (barely-lonely, Kravitz) LRC family.
  (A) finiteness bound:  every tight set has v_max <= 2n-1, the non-AP extremals hit it.
  (B) witness times: the t in (0,1) achieving gap = 1/(n+1), and their (Z/m)^* orbit.
  (C) generating rule: are non-AP tights single-element lifts v_i -> v_i + m of the AP?
  (D) Kravitz ladder: ML in {s/(ns+1)} U [1/n, inf)?  scan the gap spectrum.
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, sys
sys.path.insert(0, '04-computation')
from lrc_tight_enum_s621 import enum_tight, norm, gap_exact

def witnesses(speeds):
    """all critical t in (0,1/2] where min_i||v_i t|| == gap (the loneliness witnesses)."""
    V=[abs(v) for v in speeds]; g=gap_exact(speeds); cands=set()
    for i in range(len(V)):
        vi=V[i]
        for k in range(0,2*vi+1):
            t=Fr(2*k+1,2*vi)
            if 0<t<=Fr(1,2): cands.add(t)
        for j in range(i):
            vj=V[j]
            for d in (vi+vj,abs(vi-vj)):
                if d==0: continue
                kk=1
                while Fr(kk,d)<=Fr(1,2):
                    cands.add(Fr(kk,d)); kk+=1
    return g, sorted(t for t in cands if min(norm(v*t) for v in V)==g)

TIGHT = {n: enum_tight(n, max(2*n+4, 14), K=120) for n in range(3,8)}

print("="*78); print("(A) FINITENESS BOUND v_max <= 2n-1"); print("="*78)
ok=True
for n,sets in TIGHT.items():
    mx=max(max(s) for s in sets); bound=2*n-1
    hit=[s for s in sets if max(s)==bound]
    print(f" n={n}: v_max over all tights = {mx}  (2n-1={bound})  within bound={mx<=bound}; sets hitting 2n-1: {hit}")
    ok = ok and mx<=bound
print(" ==> bound holds for all n=3..7:", ok)

print("\n"+"="*78); print("(B) WITNESS TIMES and witness multiset (positions k where t=k/m)"); print("="*78)
for n,sets in TIGHT.items():
    m=n+1
    for s in sets:
        g,W=witnesses(list(s))
        # express witnesses with denominator m if possible
        wm=[ (w*m) for w in W]
        print(f" n={n} {str(s):24s} gap={g} witnesses(t)={[str(w) for w in W]}  (x{m}: {[str(x) for x in wm]})")

print("\n"+"="*78); print("(C) GENERATING RULE: single-element lifts v_i -> v_i + k*m of the AP that stay tight"); print("="*78)
for n in range(3,8):
    m=n+1; delta=Fr(1,m); ap=list(range(1,n+1)); gens=[]
    for i in range(n):
        for k in range(1,6):
            t=ap[:]; t[i]=ap[i]+k*m
            if reduce(gcd,t)!=1: continue
            ts=tuple(sorted(t))
            if gap_exact(t)==delta: gens.append((f"{ap[i]}->{ap[i]+k*m}", ts))
    print(f" n={n}: AP single-lifts that stay tight: {gens if gens else 'NONE (only the AP)'}")

print("\n"+"="*78); print("(D) KRAVITZ LADDER: gap spectrum vs {s/(ns+1)} and 1/n"); print("="*78)
for n in [3,4,5]:
    m=n+1; spec={}
    R=3*n+2
    for s in itertools.combinations(range(1,R+1),n):
        if reduce(gcd,s)!=1: continue
        g=gap_exact(list(s)); spec[g]=spec.get(g,0)+1
    ladder=[Fr(t,n*t+1) for t in range(1,5)]
    print(f"\n n={n} (1/n={Fr(1,n)}); ladder s/(ns+1): {[str(x) for x in ladder]}")
    below=sorted([g for g in spec if g<Fr(1,n)])
    for g in below:
        on = any(g==L for L in ladder)
        print(f"    gap={str(g):8s}={float(g):.4f}  count={spec[g]:4d}  on-ladder={on}")
    print(f"    (#distinct gaps >= 1/n: {len([g for g in spec if g>=Fr(1,n)])})")
