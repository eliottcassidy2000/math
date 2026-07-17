#!/usr/bin/env python3
"""t3_audit_hardening_kps_S128c43.py -- audit-grade referee for the hardened T3 near-pole tree:
(A) the congruence-averaging lemma: sum over one period of 1/max(1,|r(t)|) <= (2/e)(1+ln(b'/2e));
(B) the three subcases' explicit bounds vs the true near-pole sums;
(C) worst-case boundary probe (orbit hitting r=1 at t=1)."""
import sys
from math import gcd, log, pi, sin
sys.stdout.reconfigure(line_buffering=True)
LAM=1/14
def c(h): return 2*LAM if h==0 else abs(sin(2*pi*h*LAM)/(pi*h))
def inv(x,m):
    x%=m; r0,r1,s0,s1=m,x,0,1
    while r1: q=r0//r1; r0,r1=r1,r0-q*r1; s0,s1=s1,s0-q*s1
    return s0%m if r0==1 else None
print("== (A) congruence-averaging lemma referee ==")
import random
random.seed(1)
worst=0
for trial in range(200):
    bp=random.randint(3,4000)
    kap=random.randint(1,bp-1)
    e=gcd(kap,bp)
    per=bp//e
    s=0.0
    for t in range(1,per+1):
        rr=(kap*t)%bp
        rr=min(rr,bp-rr)
        if rr>0: s+=1.0/max(1,rr)
    bound=(2.0/e)*(1+log(max(1.0,bp/(2*e))))
    ratio=s/bound if bound>0 else 0
    worst=max(worst,ratio)
print("  200 random (b',kappa): max ratio (period-sum / claimed bound) = %.4f  (lemma HOLDS iff <= 1)"%worst)
print("== (B) subcase bounds vs true near-pole sums ==")
def nearpole_true_and_bound(a,b,cc,H,TMAX=8000):
    g=gcd(a,b); g0=gcd(g,cc); d=g//g0
    a2,b2=a//g,b//g; c2=cc//g0
    L=1+log(2+cc)
    true=0.0
    ia=inv(a2,b2) if b2>1 else None
    for t in range(1,TMAX):
        if ia is None: break
        rr=(-c2*t*ia)%b2
        rr=min(rr,b2-rr)
        if rr==0: continue
        for sg in (1,-1):
            h1=sg*rr if (( -c2*t - a2*sg*rr)%b2==0) else None
            if h1 is None: continue
            h2=(-c2*t-a2*h1)//b2
            if h2==0: continue
            h3=d*t
            if max(abs(h1),abs(h2),abs(h3))<=H: continue
            true+=c(h1)*c(h2)*c(h3)
    # explicit assembled bound: (2Lhat)/(pi^3 H d) + CL/(dH) + 6Lhat/(pi^3 H)  (subcases A+B+C, x2 for t<0, x2 for h2-branch)
    Lh=1+log(2+max(a,b,cc))
    bound=4*( (2*Lh)/(pi**3*H*d) + (5*Lh)/(d*H) + (6*Lh)/(pi**3*H) )
    return true,bound
for (a,b,cc) in [(307,425,541),(671,944,1413),(313,741,1531),(420,451,873)]:
    for H in (10,80,320):
        tr,bd=nearpole_true_and_bound(a,b,cc,H)
        print("  (%d,%d,%d) H=%d: true NP(h1-branch)=%.3e  explicit bound=%.3e  ratio=%.3f"%(a,b,cc,H,tr,bd,tr/bd if bd>0 else 0))
print("== (C) adversarial boundary: triple engineered so the orbit hits r=1 at t=1 ==")
# want -c' * 1 * inv(a') = 1 mod b': choose a'=1: r(t)=|(-c' t) mod b'|: pick c' = b'-1: r(1)=1
a,b,cc=1*97, 97*89, 97*(89-1)  # g=97: a'=1, b'=89, c' = 88 -> r(1)=|(-88) mod 89| = 1
g=gcd(a,b)
print("  triple (%d,%d,%d): a'=%d b'=%d"%(a,b,cc,a//g,b//g))
for H in (10,80):
    tr,bd=nearpole_true_and_bound(a,b,cc,H)
    print("   H=%d: true=%.3e bound=%.3e ratio=%.3f"%(H,tr,bd,tr/bd if bd>0 else 0))
print("DONE")
