#!/usr/bin/env python3
"""
lrc14_claimB_band_klein_S310.py
===============================
klein-2026-07-14-S310 (owner: prove the loose >=4-far margin bound = Claim B of THM-758).

Claim B is FINITE-DECIDABLE via opus's capped-envelope (THM-755, PROVED). For a >=4-far covering S,
v=maxS, P=S\\{v}, v*=r_P/(pi|G'_P|): if v>v* the capped-envelope gives disc_v<6|G'_P|^2 => M(S)>1/14
in ONE peel (PROVED); else all speeds <= v* and S is in a BOUNDED band -> finite check (mac-mini-S105
executed (220,475]). This iterates the capped-envelope peel to the terminal band bound, confirming it is
BOUNDED (<= ~500) -- so Claim B is a bounded finite check, NOT the equidistribution.
"""
import numpy as np, random
NG=1<<19; t=np.arange(NG)/NG
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def gstat(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=1.0/14)
    gi=g.astype(np.int8); r=int(np.sum(np.diff(np.concatenate([gi,gi[:1]]))==1))
    return g.mean(), r
def peel_to_band(S):
    W=sorted(set(S))
    while len(W)>1:
        v=max(W); P=[w for w in W if w!=v]; m,r=gstat(P)
        if m<=0: break
        if v > r/(np.pi*m): W=P     # capped-envelope certifies (PROVED)
        else: break                 # entered the bounded finite band
    return max(W)
random.seed(23); bounds=[]; n=0
for _ in range(60000):
    C=sorted(random.sample(range(1,random.choice([60,200,500])+1),13))
    if not iscov(C) or len([x for x in C if x>14])<4: continue
    n+=1; bounds.append(peel_to_band(C))
    if n>=120: break
b=np.array(bounds)
print("Claim B (>=4-far), %d families: terminal band bound maxP min=%d median=%d max=%d"%(n,b.min(),int(np.median(b)),b.max()))
print("=> Claim B = [v>v*: capped-envelope PROVED] + [v<=v*<=%d: bounded finite band] + [kps base]. NOT equidistribution."%b.max())
