#!/usr/bin/env python3
"""
lrc14_clustered_tail_klein_S317.py
==================================
klein-2026-07-18-S317 (owner: prove the small-killer regime and the clustered multi-killer stratum;
DFS the "169" thread for inspiration).

THM-1015: the CLUSTERED multi-killer stratum closes by INTERVAL SURVIVAL.
  V = P u K, |P| = 13-r, |K| = r in [2,6], delta = 1/14, L_P = largest interval of G_P(delta).
  If sum_i 1/k_i < L_P*(7-r)  then M(V) >= 1/14.   NO spacing hypothesis -- only size.
  Proof: |P|=13-r => M(P) >= 1/(14-r) > 1/14 by LRC(14-r) (settled) => L_P>0; the tail lemma removes at
  most 2 r delta L_P + 2 delta sum(1/k_i) from an interval of length L_P; survival needs
  sum 1/k_i < L_P (1-2 r delta)/(2 delta) = L_P (7-r) at delta = 1/14.  Ceiling r < 1/(2delta) = 7.

THE 169 CONNECTION (the DFS payoff). The telescoping balance route pays 13/14 PER killer, so j killers
cost (13/14)^j; at j=2 that is 169/196 and the chain's tightest value is (1/12)(13/14)^2 = 169/2352 =
0.071854, barely over 1/14. That 169 = 13^2 IS the two-step loss, and clustering (max barely shrinks per
step) degrades the product to ~1/2 -- which is exactly why the clustered stratum stayed open. The tail
pays NOTHING per step: killers enter only through sum 1/k_i, additively and symmetrically, so spacing is
invisible. Obstruction was multiplicative-per-step => the fix is an additive certificate.
(Other 169 sightings: MISTAKE-104's n=13 deep well {1..11,168}, 168 = 13^2-1, witness 14/169 = 14/13^2,
equioscillating, single-lift rigidity gap 1/169; THM-897's triple-beat discrepancy 22 kappa/(169 x3).)

SMALL-KILLER REGIME: killers comparable to / below the core means bounded ratio, i.e. the COMPACT case
rho < 13 -- that is HYP-7355 / THM-1014 territory, a different theorem, not a gap in this one.
"""
import numpy as np, random
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
def g(xs): return reduce(gcd,xs)
def iscov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
def Mexact(A):
    A=sorted(A); n=len(A); ds=sorted({A[i]+A[j] for i in range(n) for j in range(i,n)})
    Av=np.array(A); bn,bd=0,1
    for dd in ds:
        m=np.arange(1,dd); r=(np.outer(m,Av))%dd
        mn=np.minimum(r,dd-r).min(axis=1); val=int(mn.max())
        if val*bd>bn*dd: bn,bd=val,dd
    return Fr(bn,bd)
DLT=1.0/14
def maxL(B,dlt=DLT):
    cuts=[0.0,1.0]
    for v in B:
        iv=1.0/v; dv=dlt/v
        for k in range(0,v+1):
            for e in (k*iv-dv,k*iv+dv):
                if 0.0<=e<=1.0: cuts.append(e)
    cs=sorted(set(cuts)); best=run=0.0; prev=None
    for i in range(len(cs)-1):
        a,b=cs[i],cs[i+1]
        if b<=a: continue
        mid=0.5*(a+b); ok=True
        for v in B:
            f=(v*mid)%1.0
            if f>0.5: f=1.0-f
            if f<dlt-1e-12: ok=False; break
        if ok:
            run = run+(b-a) if (prev is not None and abs(prev-a)<1e-15) else (b-a)
            prev=b; best=max(best,run)
        else: prev=None
    return best
def criterion(P,K):
    """(star): sum 1/k_i < L_P*(7-r)  =>  M(P u K) >= 1/14."""
    r=len(K); return sum(1.0/k for k in K) < maxL(P)*(7-r)
if __name__=="__main__":
    random.seed(21)
    print("r | core | min L_P | killer threshold r/(L_P(7-r))")
    for r in range(2,7):
        Ls=[maxL(sorted(random.sample(range(1,20),13-r))) for _ in range(300)]
        Lm=min(Ls); print("%d |  %2d  | %.5f | %.1f"%(r,13-r,Lm,r/(Lm*(7-r))))
    T=Fr(1,14); n=0; bad=0; mn=Fr(1)
    for _ in range(3000):
        r=random.choice([2,3,4,5]); P=sorted(random.sample(range(1,16),13-r))
        base=random.randint(200,600)
        K=sorted(base+random.randint(-20,20) for _ in range(r))    # CLUSTERED
        V=sorted(set(P+K))
        if len(V)!=13 or g(V)!=1 or not iscov(V): continue
        n+=1; M=Mexact(V); mn=min(mn,M)
        if M<=T: bad+=1
    print("clustered census: %d families, min M = %s (%.2fx of 1/14), violations = %d"%(n,mn,float(mn)*14,bad))
