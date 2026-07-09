#!/usr/bin/env python3
"""
klein-2026-07-09-S200: does the a-priori arc-count bound c(L)<=0.37 (mac-mini-S62,
large spread) survive PRIMITIVE 7-STRUCTURED sets? (MISTAKE-128: diffs=0 mod 7 spike #arcs)

My counterexample E={0,7,14,...,82} (c=0.878) is small-spread; its dilates are
non-primitive. QUESTION: can a PRIMITIVE set with spread>=200, longest-AP<=7, and
MANY differences =0 mod 7 (gcd=1) have c=#arcs/spread > 0.37 (breaking mac-mini's
a-priori bound), or does the 7-spike vanish at large primitive spread?

Construct: a step-7 skeleton {0,7,14,...,7t} (t<=6 to keep longest-AP<=7) at scale,
PLUS off-7 points, spread up to ~500, gcd=1. Measure c and compare to 0.37.
"""
import numpy as np
from math import gcd
from functools import reduce
rng=np.random.default_rng(200)
INV7=1/7; M=6/7
def gcdl(xs): return reduce(gcd, xs)
def longest_AP(E):
    Es=sorted(set(E)); S=set(Es); b=1
    for i in range(len(Es)):
        for j in range(i+1,len(Es)):
            d=Es[j]-Es[i]; L=2; x=Es[j]+d
            while x in S: L+=1; x+=d
            b=max(b,L)
    return b
def measures(E,Nx=800000):
    E=np.array(sorted(set(E))); y=(np.arange(Nx)+0.5)/Nx
    ph=np.sort(np.mod(np.outer(y,E),1.0),axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    good=g.max(axis=1)>INV7+1e-12; gi=good.astype(int); ed=np.diff(np.concatenate([gi,gi[:1]]))
    nc=1 if gi.all() else int((ed==1).sum())
    W=np.maximum(g-INV7,0).sum(axis=1); m1,m2,m3=W.mean(),(W**2).mean(),(W**3).mean()
    return nc, good.mean(), m1/M+(m1-m2/M)**2/(m2-m3/M)

print("PRIMITIVE 7-structured, spread>=200, longest-AP<=7:  c=#arcs/spread vs 0.37 (mac-mini a-priori)")
print(f"{'spread':>7} {'#mult7':>7} {'longAP':>7} {'gcd':>4} {'#arcs':>6} {'c':>6} {'mu':>6} {'D3':>6} {'c<=0.37?':>9}")
worst=0; worstE=None
for _ in range(4000):
    # step-7 skeleton of length t (<=6) at multiplier a; then off-7 fillers; large spread
    t=int(rng.integers(3,7))           # 7-AP length t <=6 (longest-AP from this = t)
    a=int(rng.integers(1,12))          # skeleton step = 7a (still =0 mod 7)
    skel=[7*a*i for i in range(t)]
    hi=max(skel)+int(rng.integers(50,300))
    # add more mult-of-7 points (not extending the AP) + a few off-7 for primitivity
    extra7=[7*int(rng.integers(1,hi//7)) for _ in range(int(rng.integers(2,6)))]
    off7=[int(rng.integers(1,hi)) for _ in range(int(rng.integers(2,5)))]
    while any(x%7==0 for x in off7): off7=[int(rng.integers(1,hi)) for _ in off7]
    E=sorted(set([0]+skel+extra7+off7))
    if len(E)!=13: continue
    sp=max(E)
    if sp<200: continue
    if longest_AP(E)>7: continue
    if gcdl(E)!=1: continue
    nm7=sum(1 for e in E if e%7==0)
    nc,mu,d3=measures(E); c=nc/sp
    if c>worst: worst=c; worstE=(sp,nm7,longest_AP(E),nc,c,mu,d3,E)
    if c>0.37:  # report violations
        print(f"{sp:>7} {nm7:>7} {longest_AP(E):>7} {gcdl(E):>4} {nc:>6} {c:>6.3f} {mu:>6.3f} {d3:>6.3f} {'NO':>9}")
if worstE:
    sp,nm7,L,nc,c,mu,d3,E=worstE
    print(f"\nWORST: spread={sp} #mult7={nm7} longAP={L} #arcs={nc} c={c:.3f} mu={mu:.3f} D3={d3:.3f}")
    print(f"  c<=0.37? {c<=0.37}   c<mu? {c<mu}   E={E}")
print("\n=> if worst c <= 0.37 even for 7-structured primitive large-spread, mac-mini's a-priori")
print("   bound HOLDS (the 7-spike is a SMALL-spread finite-check phenomenon, covered by LEM-013).")
print("   if some c > 0.37, the a-priori #arcs bound needs a mod-7 correction term.")
