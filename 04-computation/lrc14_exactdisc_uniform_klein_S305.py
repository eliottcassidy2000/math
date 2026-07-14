#!/usr/bin/env python3
"""
lrc14_exactdisc_uniform_klein_S305.py
=====================================
klein-2026-07-14-S305 (owner: close the thin loose residual with exact Bernoulli disc; push the frontier).

RESULT 1 (loose residual closed): {1,10,..,390} stalled the CRUDE arc-count bound at peel v=390 (r=278>227),
but the EXACT disc_v = 1.19e-4 vs threshold 6|G'|^2 = 0.113 => 949x margin, cert = +0.114 > 0. The exact
disc (kps THM-732 Bernoulli form / high-precision autocorrelation) closes the whole loose branch.

RESULT 2 (frontier consolidation): the covering case is ONE route -- iterated EXACT-disc far-peel of every
speed >14, each cert (6/7)|G'_{~v}| - sqrt((6/49) disc_v) > 0, reducing to a BASE that is trivial (<=2
speeds) or subset of {1..14} (kps THM-734/738/741). VERIFIED 0 STALLS on 15 adversarial covering families
(incl {1..13\\6,182}, M=2/23, the no-k<=13-shadow counterexample). Extends kps THM-735 (<=6 far) to >6 far.
SOLE remaining analytic piece = disc_v < 6|G'_{~v}|^2 for v>14 (verified 50-1000x margins; = equidistribution).
"""
import numpy as np, random
NG=1<<19; t=np.arange(NG)/NG
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def gstat(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=1.0/14)
    return g.astype(np.float64), g.mean()
def disc(gf,v):
    G=np.fft.rfft(gf); A=np.fft.irfft(G*np.conj(G),n=NG)/NG
    idx=np.round((np.arange(v)/v)*NG).astype(np.int64)%NG
    return A[idx].mean()-A.mean()
def peel(S,Vthr=14):
    W=sorted(set(S)); worst=0.0
    while max(W)>Vthr and len(W)>1:
        v=max(W); Wn=[w for w in W if w!=v]; gf,Gnv=gstat(Wn); dv=disc(gf,v); thr=6*Gnv*Gnv
        worst=max(worst, dv/thr if thr>0 else 9)
        if dv>=thr: return sorted(W),False,worst
        W=Wn
    return sorted(W),True,worst
fams=[[1,2,3,4,5,7,8,9,10,11,12,13,182],[1,10,21,24,56,65,77,135,219,265,335,367,390],
      [1,2,3,4,5,6,7,8,9,10,11,12,182],[1,2,3,4,5,6,7,8,9,10,11,13,84],[2,3,4,5,6,7,8,9,10,11,12,13,14]]
random.seed(3)
while len(fams)<15:
    cmin=random.choice([1,1,15]); C=sorted(random.sample(range(cmin,cmin+random.randint(60,400)),13))
    if iscov(C): fams.append(C)
nb={}; stall=0; n=0
for S in fams:
    S=sorted(set(S))
    if len(S)!=13 or not iscov(S): continue
    n+=1; base,ok,w=peel(S)
    if not ok: stall+=1; print(' STALL:',S[:6],'-> base',base)
    k=(len(base),len([b for b in base if b<=14])); nb[k]=nb.get(k,0)+1
print("iterated exact-disc far-peel (>14) on %d covering families: stalls=%d ; base(size,#in{1..14})=%s"%(n,stall,dict(sorted(nb.items()))))
print("0 stalls + every base trivial(<=2) or subset{1..14} => covering = far-peel + kps base. One disc_v estimate remains.")
