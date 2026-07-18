#!/usr/bin/env python3
"""klein_synthesis_kps_S128c61.py -- kind-pasteur S128 cont.61.
(A) CROSS-VALIDATION against klein-S327 THM-1042, computed independently.
(B) SYNTHESIS: klein's (mu, N) accounting vs my largest-component accounting.
    mine:   bad <= L/7 + 2/(7k)  on one component  =>  k > 1/(3L)
    klein:  bad <= mu/7 + N/(7k) over ALL N components => k > N/(6 mu)
    Since L >= mu/N always, 1/(3L) <= N/(3 mu), so klein's form is up to 2x BETTER.
    Measure which wins in my large-killer regime.  PRINT DATA ONLY."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
    bps={F(0),F(1)}
    for p in P:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                v=F(j,p)+s
                if 0<=v<=1: bps.add(v)
    B=sorted(bps); out=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((a,b))
    mg=[]
    for a,b in out:
        if mg and mg[-1][1]==a: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return mg
print("### (A) cross-validation vs klein-S327 THM-1042 (independent computation) ###")
print("  B          my mu        klein mu    my L_max      my 1/L_max   klein 1/L_max")
KLEIN={3:(0.69048,4.2),4:(0.61905,6.2),5:(0.50476,8.8),6:(0.45714,12.0),7:(0.33469,16.3),
       8:(0.26582,22.4),9:(0.18107,31.5),10:(0.13798,46.7),11:(0.05633,77.0)}
agree=True
for k in range(3,12):
    P=list(range(1,k+1)); iv=safe_set(P)
    mu=sum(b-a for a,b in iv); L=max(b-a for a,b in iv)
    km,kl=KLEIN[k]
    ok = abs(float(mu)-km)<5e-5 and abs(float(1/L)-kl)<0.1
    agree &= ok
    print("  {1..%-2d}    %-12.5f %-11.5f %-13s %-12.1f %-12.1f %s"%(
        k,float(mu),km,str(L),float(1/L),kl,"MATCH" if ok else "*** MISMATCH"))
print("  all rows agree with klein: %s"%agree)
print()
print("### (B) klein's (mu,N) threshold  N/(6 mu)  vs mine  1/(3L) ###")
print("  setting            after removing        mu'        N     N/(6mu')   1/(3L)    winner")
def rm(iv,k):
    out=[]
    for a,b in iv:
        jlo=int(a*k)-1; jhi=int(b*k)+1; cur=a
        for j in range(jlo,jhi+1):
            x=F(j,k)-F(1,14*k); y=F(j,k)+F(1,14*k)
            if y<=a or x>=b: continue
            x=max(x,a); y=min(y,b)
            if x>cur: out.append((cur,x))
            cur=max(cur,y)
        if cur<b: out.append((cur,b))
    return out
CASES=[("r=2 worst (THM-1051)",[x for x in range(1,13) if x!=1],[873]),
       ("r=2 drop12",[x for x in range(1,12)],[873]),
       ("r=3 worst (THM-1061)",[1,2,4,5,6,7,8,9,10,11],[864,897]),
       ("r=3 alt",[2,3,4,5,6,7,8,9,10,12],[890,896]),
       ("r=4 sample",[1,2,4,5,6,7,8,10,11],[860,880,897]),
       ("r=4 sample2",[2,3,5,6,7,8,9,11,12],[870,885,899]),
       ]
for tag,P,KS in CASES:
    iv=safe_set(P)
    for k in KS: iv=rm(iv,k)
    if not iv: print("  %-19s (swallowed)"%tag); continue
    mu=sum(b-a for a,b in iv); N=len(iv); L=max(b-a for a,b in iv)
    tk=float(N/(6*mu)); tm=float(1/(3*L))
    print("  %-19s %-21s %-10.7f %-5d %-10.1f %-9.1f %s"%(
        tag,str(KS),float(mu),N,tk,tm,"klein" if tk<tm else "mine"))
print()
print("### (C) so what is the r=4 split point under the better bound? ###")
import random
random.seed(611)
C9=[sorted(c) for c in itertools.combinations(range(1,13),9)]
worst=None
for P in C9:
    iv0=safe_set(P); M=max(P); lo=13*M+1
    pool=list(range(max(lo,860),900))
    for tri in itertools.combinations(pool,3):
        iv=iv0
        for k in tri: iv=rm(iv,k)
        if not iv: continue
        mu=sum(b-a for a,b in iv); N=len(iv)
        t=float(N/(6*mu))
        if worst is None or t>worst[0]: worst=(t,tuple(P),tri)
print("  r=4 worst (mu,N) threshold over 220 cores, killer triples in [860,900): %.1f"%worst[0])
print("     at core %s killers %s"%(list(worst[1]),worst[2]))
print("DONE")
