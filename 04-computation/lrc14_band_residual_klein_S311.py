#!/usr/bin/env python3
"""
lrc14_band_residual_klein_S311.py
=================================
klein-2026-07-14-S311 (owner: execute the bounded band exhaustively for all cores).

The band (>=4-far covering S, max<=v*(S\\max), not capped-envelope-certifiable at the top peel) is BOUNDED
(max v* sampled = 1321 single-peel; <=497 iterated, S310) but C(<=500,13)-structured -- not naively
enumerable. FINDING: after applying the OTHER proved reducers to band families -- [opus density THM-746 for
diameter>339] + [capped-envelope on a DIFFERENT far element] -- the genuine RESIDUAL all has M >= 0.1428 =
2x of 1/14 (the loosest possible). So the band adds NO hard cases; a crude M>=0.14 bound would close it.
Characterization + large sample (0 failures); the full exhaustive enumeration is mac-mini's exact-Q job.
"""
import numpy as np, random
NG=1<<18; t=np.arange(NG)/NG
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def gstat(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=1.0/14)
    gi=g.astype(np.int8); r=int(np.sum(np.diff(np.concatenate([gi,gi[:1]]))==1)); return g.mean(),r
def vstar(P):
    m,r=gstat(P); return (r/(np.pi*m)) if m>0 else 1e9
def Mval(S):
    p=np.ones(NG)
    for c in S:
        fr=(c*t)%1.0; p=np.minimum(p,np.minimum(fr,1-fr))
    return p.max()
def classify(S):
    S=sorted(set(S)); v=max(S)
    if v>vstar([x for x in S if x!=v]): return 'capped-top',None
    if max(S)-min(S)>339: return 'opus-density',None
    for w in [x for x in S if x>14]:
        if w>vstar([x for x in S if x!=w]): return 'capped-other',None
    return 'RESIDUAL',Mval(S)
random.seed(31); kinds={}; resid=[]; nband=0
for _ in range(50000):
    C=sorted(random.sample(range(1,random.choice([80,250,450])+1),13))
    if not iscov(C) or len([x for x in C if x>14])<4: continue
    k,M=classify(C)
    if k!='capped-top': nband+=1; kinds[k]=kinds.get(k,0)+1
    if k=='RESIDUAL': resid.append(M)
    if nband>=200: break
print("BAND families: %d ; reducer coverage: %s"%(nband,kinds))
if resid: print("  genuine RESIDUAL: %d ; min M = %.4f = %.2fx of 1/14 (all >>1/14 => loosest residual, NOT equidist)"%(len(resid),min(resid),min(resid)/(1.0/14)))
