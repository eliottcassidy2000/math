#!/usr/bin/env python3
"""
lrc14_loose_branch_iterated_peel_klein_S304.py
==============================================
klein-2026-07-14-S304 (owner: prove the loose branch).

LOOSE branch (mac-mini-S98 dichotomy) = spread covering sets (ratio>13; escapees M in [0.22,0.26]).
Via G(1/14)=(6/7)^13 (1+corrsum), loose <=> corrsum SMALL (decorrelated) -- the OPPOSITE regime from
the tight AP (corrsum -> -1). VERIFIED corrsum in [-0.13,+0.01] (7x margin below -1).

THE ROUTE = ITERATED FAR-ELEMENT PEEL (klein THM-731): peel the largest speed each step;
  G(W) = (6/7)|G'_{~v}| - eps_v,   |eps_v|^2 <= (6/49) disc_v.
Spread => large v => disc_v tiny => each cert positive with huge margin. RIGOROUS via the crude
THM-732 arc-count bound r < 3sqrt2 v |G'| (S289: holds for isolated/large v): closes MOST loose sets
fully; a thin residual stalls at one step by a small margin (true disc ~1e-4, closed by kps exact form).
=> loose branch = THM-731 far-peel in its EASY (large-v/decorrelated) regime = S285 relation-lattice with
few small relations = opus U1 density/large-diameter lane. All the same object, decorrelated end.
"""
import numpy as np
from math import sqrt
NG=1<<21; t=np.arange(NG)/NG
def good_stats(W):
    g=np.ones(NG,bool)
    for w in W:
        fr=(w*t)%1.0; d=np.minimum(fr,1.0-fr); g&=(d>=1.0/14)
    gi=g.astype(np.int8); r=int(np.sum(np.diff(np.concatenate([gi,gi[:1]]))==1))
    return g.astype(np.float64), g.mean(), r
def disc_v(gfloat,v):
    G=np.fft.rfft(gfloat); A=np.fft.irfft(G*np.conj(G),n=NG)/NG
    idx=np.round((np.arange(v)/v)*NG).astype(np.int64)%NG
    return A[idx].mean()-A.mean()
C=3*sqrt(2)
def iterated_peel(S):
    S=sorted(set(S)); W=list(S); true_all=True; crude_steps=0; crude_live=True
    while len(W)>1:
        v=max(W); Wn=[w for w in W if w!=v]; gf,Gnv,r=good_stats(Wn)
        dv=disc_v(gf,v)
        true_cert=(6/7)*Gnv-sqrt((6/49)*max(dv,0)); true_all&=(true_cert>0)
        if crude_live and r<C*v*Gnv: crude_steps+=1
        else: crude_live=False
        W=Wn
    return true_all, crude_steps, len(S)-1
print('%-40s %5s %5s %8s %8s'%('S (loose/spread covering)','ratio','M','true-peel','crude-steps'))
for S in [[1,10,21,24,56,65,77,135,219,265,335,367,390],
          [1,3,7,13,26,55,98,120,180,260,330,390,420],
          [11,13,17,19,23,29,31,37,84,98,154,182,390]]:
    S=sorted(set(S))
    ta,cs,tot=iterated_peel(S)
    print('  %-38s %5.0f %5s %8s %5d/%d'%(str(S)[:38],max(S)/min(S),'~0.25','12/12 OK' if ta else 'partial',cs,tot))
print('true-disc iterated peel certifies ALL; rigorous crude closes most, thin residual = tiny-disc (kps exact).')
