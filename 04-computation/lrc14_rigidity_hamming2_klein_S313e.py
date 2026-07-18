#!/usr/bin/env python3
"""
lrc14_rigidity_hamming2_klein_S313e.py
======================================
klein-2026-07-17-S313e (owner: prove the inverse/rigidity theorem "M near 1/13 => A near AP").

THM-1005 (PROVED, Hamming radius 2): A = (AP\{j,k}) u {w1,w2}, w_i not in AP=[1..12]  =>  M(A) >= 2/25.
Three regimes, overlapping:
  TAIL LEMMA (THM-1004): an interval I of length L in G_S(delta) survives r new speeds once
      2 r delta L + 2 delta sum(1/w_i) < L.
  R1 both large : r=2, delta=2/25 -> sum 1/w_i < 17L/4; holds if both w > 8/(17 L_{j,k}).
                  W_joint = max over 66 pairs = 8000/119 = 67.23
  R2 mixed      : fix w1<=68, B'=B u {w1} (11 speeds, M>=1/12>delta), r=1 -> w2 > 4/(21 L').
                  W2 = 11400/217 = 52.53 (worst j,k=4,6, w1=38); no degenerate L'=0.
  R3 both small : exact check, 13<=w1<w2<=68, all 66 (j,k) = 101,640 families, ZERO violations, min=2/25.
Since w>=13 and 68 > max(W_joint,W2), the regimes exhaust all (w1,w2).
Extremals (sharp, and NOT unique): {1..11,24} (radius 1) and {1,2,3,5,7,8,9,10,11,12,17,19} (radius 2).
"""
import numpy as np, itertools
from fractions import Fraction as Fr
AP=list(range(1,13)); delta=Fr(2,25); num,den=2,25
def M_ge_delta(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A)
    for d in ds:
        m=np.arange(1,d); r=(np.outer(m,Av))%d
        mn=np.minimum(r,d-r).min(axis=1)
        if (mn*den >= num*d).any(): return True
    return False
def maxL(B,d=delta):
    cuts={Fr(0),Fr(1)}
    for v in B:
        for k in range(0,v+1):
            for e in (Fr(k,v)-d/v, Fr(k,v)+d/v):
                if 0<=e<=1: cuts.add(e)
    cs=sorted(cuts); segs=[]
    for i in range(len(cs)-1):
        a,b=cs[i],cs[i+1]
        if b>a and all(min((v*((a+b)/2))%1,1-(v*((a+b)/2))%1)>=d for v in B): segs.append([a,b])
    m=[]
    for a,b in segs:
        if m and m[-1][1]==a: m[-1][1]=b
        else: m.append([a,b])
    return max((b-a for a,b in m), default=Fr(0))
if __name__=="__main__":
    Wj=max(Fr(8,17)/maxL([v for v in AP if v not in (j,k)]) for j,k in itertools.combinations(AP,2))
    print("R1 W_joint =",Wj,"=%.2f"%float(Wj))
    W2=Fr(0)
    for j,k in itertools.combinations(AP,2):
        B=[v for v in AP if v not in (j,k)]
        for w1 in range(13,69): W2=max(W2,Fr(4,21)/maxL(B+[w1]))
    print("R2 W2      =",W2,"=%.2f"%float(W2))
    bad=0; n=0
    for j,k in itertools.combinations(AP,2):
        B=[v for v in AP if v not in (j,k)]
        for w1,w2 in itertools.combinations(range(13,69),2):
            n+=1
            if not M_ge_delta(B+[w1,w2]): bad+=1
    print("R3 finite  : %d families, violations=%d"%(n,bad))
    print("=> THM-1005 PROVED" if bad==0 and float(Wj)<=68 and float(W2)<=68 else "=> NOT closed")
