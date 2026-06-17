#!/usr/bin/env python3
"""
lrc14_compactness_probe_kps  (kind-pasteur, PROVE side)

COMPACTNESS LEMMA probe: is it true that any primitive 13-set with small L
must be a BOUNDED-lcm perturbation of a tight config?

Strategy: random + structured search for 13-sets (distinct positive ints, gcd 1)
with small L, recording for each:
  - L (exact)
  - max element, lcm, "core" = the part agreeing with {1..13}
  - hamming distance from {1..13} (how many elements differ)
We collect ALL found with L < threshold and inspect: do they all have bounded
max-element / are they all <=1-element or 2-element perturbations of a tight set?

This is EMPIRICAL evidence for/against the compactness lemma. We CANNOT prove
finiteness by search, but we can (a) confirm the small-L locus is sparse and
clustered near tight configs in the searched range, and (b) find the SECOND-
smallest L values to see the spectrum above 1/1260.
"""
from fractions import Fraction as F
import math, random

def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=max((14*k-1)*inv,0.0); hi=min((14*k+1)*inv,1.0)
            if lo<hi: arcs.append((lo,hi))
    arcs.sort(); tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=hi if hi>ch else ch
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot
def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(x,y) for x,y in out if y>x]
def L_exact(S):
    arcs=[]
    for v in set(S): arcs+=danger(v)
    arcs.sort(); t=F(0); cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=ch if ch>hi else hi
        else: t+=ch-cl; cl,ch=lo,hi
    return F(1)-(t+ch-cl)
def lcm(S):
    l=1
    for v in S: l=l*v//math.gcd(l,v)
    return l
def ham_from_AP(S):
    return 13-len(set(S)&set(range(1,14)))

TIGHTS=[frozenset(range(1,14)), frozenset(list(range(1,12))+[13,24])]
THRESH=1.0/100.0   # 0.01
found={}  # config-> L_exact

# 1) all single perturbations of both tight configs, w<=2000
print("[1] single perturbations of tight configs (w<=2000)",flush=True)
for T in TIGHTS:
    base=sorted(T)
    for e in base:
        rest=[x for x in base if x!=e]
        for w in range(2,2001):
            if w in rest: continue
            S=tuple(sorted(rest+[w]))
            if len(S)!=13: continue
            if L_float(S)<THRESH:
                found.setdefault(S,None)

# 2) all 2-element perturbations of AP, w<=80 (float screen)
print("[2] 2-element perturbations of AP (w<=80)",flush=True)
from itertools import combinations
base=list(range(1,14))
for e1,e2 in combinations(base,2):
    rest=[x for x in base if x not in (e1,e2)]
    for w1,w2 in combinations(range(2,81),2):
        if w1 in rest or w2 in rest: continue
        S=tuple(sorted(rest+[w1,w2]))
        if len(S)!=13: continue
        if L_float(S)<THRESH:
            found.setdefault(S,None)

# 3) random 13-sets with bounded max element (broad net)
print("[3] random 13-sets max<=60, 300k samples",flush=True)
random.seed(1)
for _ in range(300000):
    S=tuple(sorted(random.sample(range(1,61),13)))
    if L_float(S)<THRESH:
        found.setdefault(S,None)

print(f"\nfound {len(found)} configs with L_float < {THRESH}. Exact-evaluating...",flush=True)
rows=[]
for S in found:
    Le=L_exact(list(S))
    if Le<F(1,100):
        rows.append((Le,S))
rows.sort(key=lambda r:(r[0],))
print(f"\nConfigs with EXACT L < 1/100, sorted (showing up to 40):")
print(f"  {'L':>14s}  {'float':>10s} {'maxel':>5s} {'lcm':>8s} {'hAP':>4s}  config")
for Le,S in rows[:40]:
    print(f"  {str(Le):>14s} {float(Le):>10.3e} {max(S):>5d} {lcm(S):>8d} {ham_from_AP(S):>4d}  {list(S)}")

print(f"\nTotal with L<1/100: {len(rows)}")
# spectrum of distinct L values
vals=sorted(set(Le for Le,_ in rows))[:15]
print("Smallest distinct L values found:")
for v in vals:
    print(f"   {v} = {float(v):.6e}")
# max element among small-L configs
if rows:
    print(f"\nMax element seen among L<1/100 configs: {max(max(S) for _,S in rows)}")
    print(f"Max lcm seen: {max(lcm(S) for _,S in rows)}")
    print(f"Max hamming-from-AP: {max(ham_from_AP(S) for _,S in rows)}")
