#!/usr/bin/env python3
"""
lrc14_scaling_proof_kps-S6-wf.py
Bulletproof verification of the SCALING INVARIANCE THEOREM (the one rigorous
reduction this angle establishes), plus the reflection and a clean restatement.

THEOREM (Scaling).  For any finite E subset Z and any positive integer c,
   meas(N(cE)) = meas(N(E)),   where N(E)={x in [0,1): maxgap{frac(ex):e in E}<=1/7}.
   More generally mu_theta(cE)=mu_theta(E) for any theta.
PROOF.  frac(c e x)=frac(e (cx)) since e in Z.  The map phi(x)=cx mod 1 is a
   measure-preserving c-to-1 endomorphism of [0,1).  The gap-configuration of
   {frac(c e_i x)} equals that of {frac(e_i y)} at y=phi(x).  Hence the indicator
   1[maxgap<=1/7] of cE at x equals that of E at phi(x); integrating and using
   measure-preservation gives equality.  QED.

COROLLARY.  WLOG gcd(E)=1 (divide by gcd).  So the spread bound need only be
   proved for PRIMITIVE E.  (But note: primitivity does NOT bound spread; the
   single far point [0..6, 14] is primitive with arbitrarily large spread and
   net>0.  So a finite reduction does NOT follow from scaling alone.)

This script verifies the theorem on a stress battery and prints the exact
consecutive values mu_1/7(consec_k) used downstream.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(2024)
ONE7=F(1,7)

def mu_theta(E, theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[(E[order[t]]*mid).__floor__() for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; k1=ks[t]; k2=ks[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s);  subs.append((lo,b)) if lo<b else None
            else:
                hi=min(b,(theta-c)/s);  subs.append((a,hi)) if a<hi else None
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

if __name__=="__main__":
    print("[SCALING THEOREM] mu_theta(cE)==mu_theta(E) stress battery:")
    bad=0; tot=0
    for _ in range(120):
        k=random.choice([6,7,8,9,10])
        E=sorted(set([0]+random.sample(range(1,18),k-1)))
        if len(E)<k: continue
        for theta in [F(1,7),F(2,7),F(1,5),F(3,11)]:
            base=mu_theta(E,theta)
            for c in (2,3,4,5,6,7):
                tot+=1
                if mu_theta([c*e for e in E],theta)!=base:
                    bad+=1; print(f"   FAIL E={E} c={c} theta={theta}")
    print(f"   {tot} checks across theta in {{1/7,2/7,1/5,3/11}}, c in 2..7: violations={bad}")
    print(f"   => SCALING INVARIANCE {'VERIFIED (consistent with the proof)' if bad==0 else 'FALSIFIED'}")

    print("\n[CONSEC VALUES] mu_1/7(consecutive_k), the downstream floor inputs:")
    for k in range(7,14):
        v=mu_theta(list(range(k)),ONE7)
        print(f"   k={k:2d}: mu_1/7(consec)={v}={float(v):.6f}  net=1-mu={1-v}={float(1-v):.6f}")
