#!/usr/bin/env python3
"""
lrc14_general_arcwidth_criterion — mac-mini-2026-06-17-S5

THE GENERALIZED ARC-WIDTH LEMMA (any runner, not just multiples of 14):
For ANY speed v, the danger set {tau: ||v tau||<1/14} is v TEETH of full width 1/(7v),
centers k/v spaced 1/v = 7*(width) apart (gaps 6/(7v)). So if the level-1/14 safe set of
the OTHER runners S\{v} has a widest arc W(S\{v}) > 1/(7v), that arc cannot fit in a tooth
=> it contains a v-safe point => M(S) >= 1/14.

CRITERION C(S):  EXISTS v in S with  W(S\{v}) > 1/(7v).   C(S) => M(S) >= 1/14.

THE PROOF QUESTION: does EVERY covering 13-set satisfy C(S)? If yes, LRC(14) is PROVED
(for the only hard case). Test over structured + random covering sets, including configs
with TWO large runners. For each report whether some v works, which v, and the margin.
"""
from fractions import Fraction as F
from itertools import combinations
import random

C=F(1,14)
def darcs(v,c=C):
    hw=F(c,v); return [(F(k,v)-hw,F(k,v)+hw) for k in range(v)]
def wrapU(iv):
    o=[]
    for lo,hi in iv:
        s=lo-(lo%1); a=lo-s;b=hi-s
        if b<=1:o.append((a,b))
        else:o.append((a,F(1)));o.append((F(0),b-1))
    o=sorted(o);r=[];cl,ch=o[0]
    for lo,hi in o[1:]:
        if lo<=ch: ch=ch if ch>hi else hi
        else:r.append((cl,ch));cl,ch=lo,hi
    r.append((cl,ch));return r
def Wsafe(A,c=C):
    """widest arc of the level-c safe set of A (complement of danger teeth)."""
    dz=[]
    for v in set(A): dz+=darcs(v,c)
    if not dz: return F(1)
    dz=wrapU(dz)
    # if danger covers everything, no safe arc
    if sum(hi-lo for lo,hi in dz)>=1:
        # still may have measure-zero/positive gaps; compute gaps
        pass
    best=F(0)
    for i in range(len(dz)):
        hi=dz[i][1]; lo=dz[(i+1)%len(dz)][0]+(1 if i==len(dz)-1 else 0)
        if lo-hi>best: best=lo-hi
    return best
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def criterion(S):
    """return (holds, best_v, best_margin) for C(S): exists v with W(S\{v})>1/(7v)."""
    best=None
    for v in sorted(set(S)):
        A=[u for u in S if u!=v]
        W=Wsafe(A); thr=F(1,7*v); margin=W-thr
        if best is None or margin>best[2]: best=(margin>0,v,margin,W,thr)
        if margin>0: return (True,v,margin,W,thr)
    return (False,)+best[1:]

# exact M for cross-check
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); Cc=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cc.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cc.add(F(k,d)); k+=1
    Cc.add(F(1,2)); return Cc
def M(S): return max(g(S,t) for t in cand(S))

print("="*78)
print("CRITERION C(S): exists v with W(S\\{v}) > 1/(7v)  =>  M(S)>=1/14")
print("="*78)
named=[
 ("min {1..11,13,84}",[1,2,3,4,5,6,7,8,9,10,11,13,84]),
 ("drop6 {1..5,7..13,84}",[1,2,3,4,5,7,8,9,10,11,12,13,84]),
 ("drop6 u98",[1,2,3,4,5,7,8,9,10,11,12,13,98]),
 ("two-large {1..10,11,84,182}? ",None),
]
for name,S in named:
    if S is None: continue
    if len(set(S))!=13 or not covering(S):
        print(f"  {name}: (not a covering 13-set, skip)"); continue
    h=criterion(S); Mv=M(S)
    print(f"  {name:26s}: C holds={h[0]} via v={h[1]} (W={float(h[3]):.5f} > 1/(7v)={float(h[4]):.5f}, margin={float(h[2]):+.5f})  M(S)={Mv}={float(Mv):.4f}")

print("\n"+"="*78)
print("STRESS TEST: does C(S) hold for ALL covering 13-sets?  (incl. multiple large runners)")
print("="*78)
random.seed(5); ntot=0; cfail=0; mbreak=0; worstmargin=(F(99),None)
fails=[]
def rand_covering():
    # base small set + 1 or 2 large covering elements
    drop=random.choice([1,2,3,4,5,6,7,12])
    base=[v for v in range(1,14) if v!=drop]
    if random.random()<0.5:
        S=base+[84*random.randint(1,12)]
    else:
        # two large: remove another small, add two large covering mults
        d2=random.choice([x for x in base if x in (7,8,9,10,11,12,13)])
        base2=[v for v in base if v!=d2]
        S=base2+[84*random.randint(1,6), random.choice([182,168,252,364,154,286])*random.randint(1,3)]
    return sorted(set(S))
while ntot<3000:
    S=rand_covering()
    if len(S)!=13 or not covering(S): continue
    ntot+=1
    h=criterion(S);
    if not h[0]:
        cfail+=1; fails.append(S)
        if M(S)<C: mbreak+=1
    else:
        if h[2]<worstmargin[0]: worstmargin=(h[2],S,h[1])
print(f"  covering 13-sets tested: {ntot}")
print(f"  C(S) FAILED (no v with W>1/(7v)): {cfail}")
print(f"  of those, actual LRC break M(S)<1/14: {mbreak}")
if cfail==0:
    print("  => C(S) HOLDS UNIVERSALLY on the sample: the generalized arc-width criterion")
    print("     PROVES M(S)>=1/14 for every covering 13-set tested. (Proof target: show C always holds.)")
    print(f"  tightest successful margin W-1/(7v) = {float(worstmargin[0]):.6f} at v={worstmargin[2]} S={worstmargin[1]}")
else:
    print(f"  C failed on {cfail} sets (e.g. {fails[:2]}) — criterion not universal; need a fallback there.")
