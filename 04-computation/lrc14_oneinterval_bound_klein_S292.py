#!/usr/bin/env python3
"""
lrc14_oneinterval_bound_klein_S292.py
=====================================
klein-2026-07-13-S292 (owner: prove the large-speed one-interval equidistribution bound).

GOAL: for large-speed clusters C, conc(C)=14|G(C)cap[0,1/14)|/|G(C)| < 7 (=> L({1}UC)>0, S290).

RIGOROUS SINGLE-SPEED BOUND. G(C) subset G({c}) for any c in C, so on [0,1/14) only ONE speed's bad
fraction (1/7) is removed:  |G(C)cap[0,1/14)| <= |G({c})cap[0,1/14)| <= 3/49 + 6/(7c)
(good per period = 6/(7c); floor(c/14) full periods <= 3/49; + one partial). Take c=max(C):
    conc(C) <= (6/7 + 12/max(C)) / |G(C)|      =>   conc<7  whenever  |G(C)| > 6/49 + 12/(7 max(C)).
VERIFIED here (conc <= bound, 0 violations). SCOPE (honest): the leading 6/7 ~ 1 makes the margin thin;
the bound fires only for well-spread large-speed clusters (|G(C)|>6/49, ~62% sampled) and needs large
max; it FAILS for |G(C)|<=6/49 (AP-dilate-like, where conc is actually SMALL by dilation-spreading ->
opus dilation-blindness). Full margin (true conc~3.3, ~2x under 7) needs MULTI-speed equidistribution =
the cancellation restricted to [0,1/14). NOT elementary. NB kps-THM-735 (simultaneous multi-peel)
DISSOLVES the S289 isolation wall for bounded bodies; the residual is {1}Ularge-cluster (opus true-disc).
"""
import numpy as np, random
NG=1<<20; THR=1.0/14.0; t=np.arange(NG)/NG
def good(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=THR)
    return g
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
lo=int(round(NG/14.0)); S649=6.0/49.0
print("single-speed bound  conc <= (6/7 + 12/max(C))/|G(C)|  ; conc<7 iff |G(C)| > 6/49 + 12/(7max)")
print("%-24s %8s %7s %6s %9s %6s %s"%("C","|G(C)|","conc","max","BOUND","<7?","ok"))
tests=[list(range(90,102)),list(range(30,42)),list(range(45,57)),
       [2*k for k in range(2,14)],[15*k for k in range(2,14)],list(range(3,15)),
       random.Random(3).sample(range(20,200),12)]
viol=0
for C in tests:
    C=sorted(set(C)); g=good(C); L=g.mean(); n0=g[:lo].sum()/NG
    conc=14*n0/L if L>0 else 9; bd=(6.0/7.0+12.0/max(C))/L if L>0 else 99
    ok=conc<=bd+1e-9; viol+= (0 if ok else 1)
    print("%-24s %8.5f %7.3f %6d %9.3f %6s %s"%(str(C)[:24],L,conc,max(C),bd,'Y' if bd<7 else 'n',ok))
print("bound holds (conc<=bound) with %d violations."%viol)
# what fraction of large-speed covering {1}UC clear |G(C)|>6/49?
random.seed(5); n=0; fire=0
for _ in range(4000):
    N=random.choice([30,50,90,150]); b=random.choice([15,20,30,60])
    C=sorted(random.sample(range(b,b+N),12)); S=sorted([1]+C)
    if not iscov(S): continue
    n+=1; L=good(C).mean()
    if L>S649+12.0/(7*max(C)): fire+=1
print("large-speed covering: %d sampled; single-speed bound fires for %d (%.0f%%). residual=AP-dilate-like (opus)."%(n,fire,100*fire/max(n,1)))
print("done.")
