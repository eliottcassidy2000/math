#!/usr/bin/env python3
"""
lrc14_floor_sequences — mac-mini-2026-06-18-S3

The integer sequences of the LRC(14) S3 floor. Hunt closed forms / recurrences / OEIS.
  mu(E) = meas{ x in [0,1) : {frac(e x): e in E} have circular max-gap > 2/7 },  0 in E.
Verified exact good-set engine (piecewise-linear, validated vs mu_3=1, mu_4=19/21, mu_5=9/14,
mu_13=829/4620 from THM-527).

Sequences:
  (A) mu_consec(k) = mu({0,1,...,k-1})  -- the AP/Steinhaus floor (consecutive cluster).
  (B) F(k) = P(k iid uniform pts have max-gap > 2/7)  = 1 - sum_j (-1)^j C(k,j)(1-2j/7)_+^{k-1}
      (the equidistribution CEILING; denom 7^{k-1}).
  (C) mu_min(k) true minima (k<=7 exact from agents): 1, 19/21, 9/14, 13/35, ...
Goal: a closed form for mu_consec(k) (Farey-structured) that proves a positive floor and
maybe a clean recurrence -> simplify the spread-bound lemma.
"""
from fractions import Fraction as F
from math import comb, gcd
THR=F(2,7)
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def good_set_exact(E):
    es=sorted(set(E)); bps={F(0),F(1)}
    for e in es:
        for n in range(1,e): bps.add(F(n,e))
    for i in range(len(es)):
        for j in range(i+1,len(es)):
            d=es[j]-es[i]
            for n in range(1,d): bps.add(F(n,d))
    bp=sorted(bps); good=[]
    for t in range(len(bp)-1):
        lo,hi=bp[t],bp[t+1]; mid=(lo+hi)/2; vm={}
        for e in es:
            v=(e*mid)%1
            if v not in vm: vm[v]=(e,(e*mid)-v)
        keys=sorted(vm); se=[vm[v][0] for v in keys]; sc=[vm[v][1] for v in keys]; m=len(keys)
        if m==1: good.append((lo,hi)); continue
        sub=[]
        for r in range(m):
            if r<m-1: al=se[r+1]-se[r]; be=-(sc[r+1]-sc[r])
            else: al=se[0]-se[m-1]; be=F(1)-(sc[0]-sc[m-1])
            if al==0:
                if be>THR: sub.append((lo,hi))
            else:
                xb=(THR-be)/al
                if al>0:
                    a=max(lo,xb)
                    if a<hi: sub.append((a,hi))
                else:
                    b=min(hi,xb)
                    if lo<b: sub.append((lo,b))
        for a,b in merge(sub):
            if a<b: good.append((a,b))
    return merge(good)
def mu(E): return sum((b-a for a,b in good_set_exact(E)), F(0))
def Fiid(k):
    s=F(0); j=0
    while 1-j*THR>0:
        s+=(-1)**j*comb(k,j)*(1-j*THR)**(k-1); j+=1
    return 1-s

print("="*72)
print("(A) mu_consec(k) = mu({0..k-1})  [AP/Steinhaus floor]")
print("="*72)
muc={}
for k in range(3,23):
    m=mu(list(range(k))); muc[k]=m
    print(f"  k={k:2d}: mu={str(m):>14}  = {float(m):.6f}   num={m.numerator} den={m.denominator}")
print("\n  numerators:", [muc[k].numerator for k in range(3,23)])
print("  denominators:", [muc[k].denominator for k in range(3,23)])
# denominator factorizations
def fac(n):
    n=abs(n); f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
print("  denom factorizations:")
for k in range(3,16): print(f"    k={k:2d}: {muc[k].denominator} = {fac(muc[k].denominator)}")

print("\n"+"="*72)
print("(B) F(k) iid CEILING  (denom 7^{k-1})")
print("="*72)
for k in range(3,16):
    f=Fiid(k); print(f"  k={k:2d}: F={str(f):>22} = {float(f):.6f}  num={f.numerator}")
print("  iid numerators (over 7^{k-1}):", [(Fiid(k)*7**(k-1)).numerator if False else (1-Fiid(k))  for k in range(3,8)])

print("\n"+"="*72)
print("(C) differences / recurrence probes for mu_consec(k)")
print("="*72)
# successive ratios and differences
for k in range(3,21):
    if k+1 in muc:
        d=muc[k]-muc[k+1]
        print(f"  mu({k})-mu({k+1}) = {str(d):>16} = {float(d):.6f}")
print("\n  14*mu_consec(k) (clearing the /7 in the threshold):")
for k in range(3,16):
    print(f"    k={k:2d}: 14*mu = {14*muc[k]} ")
