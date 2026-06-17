"""
DISPROVE-side CRUX ATTACK: can meas(G_C) (lonely measure of a 12-speed set C) -> 0?

The PROVE side's open gap: a UNIFORM positive lower bound on meas(G_C) over ALL
12-subsets C that appear in a near-tight 13-set. If we can drive meas(G_C) -> 0 with
lcm(C) -> inf, then adding one speed could give L(C u {w}) -> 0 and inf L = 0,
killing the singular-series C'(14) route.

This is THE crux. We search for 12-sets C of distinct positive integers with the
SMALLEST possible lonely measure meas(G_C) = L(C). The danger arcs of 12 speeds have
total measure 12/7 ~ 1.714 > 1, so they CAN cover [0,1] (meas=0) -- the question is
whether they cover with lcm -> inf (a genuinely new tight 12-locus) or only near {1..12}.

Strategy:
 (1) Exhaustive-ish: which 12-sets are TIGHT (meas=0)? Census small windows.
 (2) Minimize meas(G_C) over random + structured 12-sets with growing max entry;
     track whether the min meas trends to 0 as we allow larger entries / lcm.
 (3) The KEY adversarial family: 13-set = C u {w}. If C is a NEW tight 12-set with
     large lcm, then a minimal perturbation gives small L with large lcm.
     We hunt tight 12-sets beyond {1..12} and its relatives.
"""
from fractions import Fraction as Fr
from math import gcd, lcm
from functools import reduce
from itertools import combinations
import random

def danger_arcs(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A
_cache={}
def darcs(v):
    r=_cache.get(v)
    if r is None: r=danger_arcs(v); _cache[v]=r
    return r
def L_exact(S):
    A=[]
    for v in S: A.extend(darcs(v))
    A.sort(key=lambda t:(t[0],t[1]))
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if b<=a: continue
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

print("=== (1) TIGHT 12-sets census in small windows ===")
# A 12-set C is tight iff its danger arcs cover [0,1]. Census subsets of [1,N].
for N in [13,14,15,16]:
    tights=[]
    for C in combinations(range(1,N+1),12):
        if reduce(gcd,C)!=1: continue
        if L_exact(C)==0:
            tights.append(C)
    print(f"   primitive tight 12-sets in [1,{N}]: {len(tights)}")
    for C in tights[:20]:
        print(f"      {C}  lcm={lcm(*C)}")

print()
print("=== (2) Minimize meas(G_C) over 12-sets, growing max entry ===")
random.seed(7)
for cap in [13,16,20,30,50,100]:
    best=None; bestC=None
    # include the obvious {1..12}
    trials=[tuple(range(1,13))]
    for _ in range(40000):
        trials.append(tuple(sorted(random.sample(range(1,cap+1),12))))
    seen=set()
    for C in trials:
        if C in seen: continue
        seen.add(C)
        if reduce(gcd,C)!=1: continue
        L=L_exact(C)
        if best is None or L<best:
            best=L; bestC=C
    print(f"   max-entry<={cap:3d}: min meas(G_C) = {float(best):.8f} (={best}) at {bestC} lcm={lcm(*bestC)}")

print()
print("=== (3) Hunt NEW tight 12-sets with large lcm (1- and 2-perturbations of {1..12}) ===")
base=list(range(1,13))
found=set()
# 1-replacement
for i in range(12):
    for w in range(13,400):
        if w in base: continue
        C=tuple(sorted(base[:i]+base[i+1:]+[w]))
        if reduce(gcd,C)!=1 and len(set(C))==12: continue
        if len(set(C))!=12: continue
        if L_exact(C)==0:
            found.add(C)
# 2-replacement (windowed)
for i,j in combinations(range(12),2):
    rem=[base[k] for k in range(12) if k not in (i,j)]
    for w1 in range(13,80):
        for w2 in range(w1+1,80):
            C=tuple(sorted(rem+[w1,w2]))
            if len(set(C))!=12: continue
            if L_exact(C)==0:
                found.add(C)
print(f"   new tight 12-sets found (excluding pure {{1..12}}): {len(found - {tuple(range(1,13))})}")
for C in sorted(found):
    print(f"      {C}  lcm={lcm(*C)}  {'(=base)' if C==tuple(range(1,13)) else ''}")
print("   => if the ONLY tight 12-sets are {1..12} and near-relatives (bounded lcm),")
print("      then meas(G_C) is bounded below on every OTHER 12-core => uniform floor plausible.")
