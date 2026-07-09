#!/usr/bin/env python3
"""
klein-2026-07-09-S196: the DISSOCIATED branch (longest-AP <= k-6). Is j* small,
and is there an ELEMENTARY mechanism (to complement LEM-012's L>=k-5 branch)?

Construct HARD (V/7-dense at j=1) sets with SHORT longest-AP and measure j*.
Test candidate elementary bounds. If j* <= small always, look for the reason.
"""
import numpy as np
from math import gcd
rng=np.random.default_rng(196963)
INV7=1/7
def maxgap(pts,V):
    p=np.sort(np.array(pts)%V); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()/V
def is_good(E,j,V): return maxgap([(e*j)%V for e in E],V)>INV7+1e-12
def jstar(E,V,Jmax=200):
    for j in range(1,min(Jmax,V)):
        if is_good(E,j,V): return j
    return None
def longest_AP(E):
    Es=sorted(set(E)); Eset=set(Es); best=1
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            best=max(best,L)
    return best

# find HARD sets with SHORT longest-AP (dissociated)
print("HARD dissociated sets (longest-AP <= k-6), j* distribution:")
print(f"{'k':>3} {'V':>5} {'#found':>7} {'maxL':>5} {'max j*':>7} {'j* hist (1..)':>20}")
for k in (8,11,13):
    for V in (91, 200, 500):
        found=[]; tries=0
        while len(found)<300 and tries<400000:
            tries+=1
            rest=rng.choice(np.arange(1,V),k-1,replace=False)
            E=tuple(sorted([0]+[int(x) for x in rest]))
            if max(E)<6*V/7: continue
            if is_good(E,1,V): continue      # hard
            if longest_AP(E)>k-6: continue   # dissociated
            found.append(E)
        if not found:
            print(f"{k:>3} {V:>5} {0:>7}  (none)"); continue
        js=[jstar(E,V) for E in found]
        mxj=max(j for j in js if j)
        maxL=max(longest_AP(E) for E in found)
        from collections import Counter
        hist=dict(sorted(Counter(js).items())[:6])
        print(f"{k:>3} {V:>5} {len(found):>7} {maxL:>5} {mxj:>7} {str(hist):>20}")

# candidate elementary mechanism: does j=2 or j=3 always work for hard dissociated?
print("\nCandidate: for hard dissociated, is j* <= 7 always? (and what's the worst)")
worst=0; worstE=None
for k in (8,11,13):
    for V in (91,200,500,1001):
        tries=0; cnt=0
        while cnt<400 and tries<500000:
            tries+=1
            rest=rng.choice(np.arange(1,V),k-1,replace=False)
            E=tuple(sorted([0]+[int(x) for x in rest]))
            if max(E)<6*V/7 or is_good(E,1,V) or longest_AP(E)>k-6: continue
            cnt+=1
            js=jstar(E,V)
            if js and js>worst: worst=js; worstE=(k,V,E)
print(f"  worst j* over hard dissociated = {worst}  at k,V={worstE[:2] if worstE else '-'}")
print(f"  worst E = {worstE[2][:8] if worstE else '-'}...")
