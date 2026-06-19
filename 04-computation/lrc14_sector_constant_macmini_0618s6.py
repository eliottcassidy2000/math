#!/usr/bin/env python3
"""
lrc14_sector_constant_macmini_0618s6.py  (mac-mini-2026-06-18-S6)

Pin the certificate constant for meas(S7(E)) <= cap_k.
Theory: corr(E) = meas(S7(E)) - M7(k) is the offset relation-lattice tail; the dominant
(support-3) part is bounded by sum over primitive triples of (sector-Fourier const)/height
~ C * W(E), W(E) = sum_{triples} 1/height (the THM-503 7-vanishing kills 7|n modes).

TESTS over a broad bank (consec, perforated, embedded-AP, dissociated, random, k=8..11):
 (1) the max ratio C* = max corr(E)/W(E) (so corr(E) <= C* * W(E) empirically).
 (2) does consec MAXIMIZE W(E)?  (then corr<=C*W<=C*W(consec))
 (3) certificate-only closure: is C* * W(consec_k) < margin_k = cap_k - M7(k)? (=> no finite
     residual needed) OR is consec at/over the boundary (=> finite low-height check needed)?
 (4) the residual {W(E) > margin_k/C*}: how many shapes, are they all AP-rich (=> finite mod scale)?
"""
import sys, itertools, random
from math import comb, gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(6186)

def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(int(((e*xm)%1)*7) for e in E)
        if len(secs)==7: total+=x1-x0
    return total
def M7(k): return sum(F((-1)**t*comb(6,t))*F(7-t,7)**(k-1) for t in range(7))
H=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H/u)%1; b=(c+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgmerge(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    if not P: return F(1)
    dz=mgmerge([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
import functools
@functools.lru_cache(None)
def cap(k):
    psz=13-k
    if psz==0: return F(1)
    return min(measGP(P) for P in itertools.combinations(range(1,14),psz))
def W(E):
    E=sorted(set(E)); w=0.0
    for i,j,l in itertools.combinations(range(len(E)),3):
        a,b,c=E[j]-E[l],E[l]-E[i],E[i]-E[j]
        g=gcd(gcd(abs(a),abs(b)),abs(c)) or 1
        h=abs((a//g)*(b//g)*(c//g))
        if h>0: w+=1.0/h
    return w

print("="*76)
print("(2,3) consec W and certificate-only closure check, k=8..11")
print("="*76)
for k in range(8,12):
    Wc=W(list(range(k))); m=float(M7(k)); c=float(cap(k)); margin=c-m
    # find max W over shapes spread<=k+2 (proxy for 'consec maximizes W')
    mxW=(Wc,'consec')
    for rest in itertools.combinations(range(1,k+3),k-1):
        E=[0]+list(rest); ww=W(E)
        if ww>mxW[0]: mxW=(ww,E)
    print(f"  k={k}: W(consec)={Wc:.3f}  maxW(spread<=k+2)={mxW[0]:.3f} "
          f"{'(consec)' if mxW[1]=='consec' else '(NOT consec: '+str(mxW[1])+')'}  margin_k={margin:.4f}")

print("\n"+"="*76)
print("(1) max ratio C* = corr/W over a broad bank (k=8); corr(E)<=C*W?")
print("="*76)
M78=float(M7(8)); bank=[]
# structured
for rest in itertools.combinations(range(1,13),7):
    bank.append([0]+list(rest))
# random spread
for _ in range(200):
    sp=random.randint(8,60); bank.append(sorted(set([0,sp]+random.sample(range(1,sp),6))))
maxC=0.0; argC=None; corr_consec=measS7(list(range(8)))-M78
for E in bank:
    if len(set(E))!=8: continue
    corr=float(measS7(E))-M78; w=W(E)
    if w>1e-9:
        r=corr/w
        if r>maxC: maxC=r; argC=E
print(f"  max corr/W = {maxC:.5f} at E={argC}")
print(f"  corr(consec_8)={float(corr_consec):.5f}, W(consec_8)={W(list(range(8))):.3f}, "
      f"ratio={float(corr_consec)/W(list(range(8))):.5f}")
margin8=float(cap(8))-M78
print(f"  C* * W(consec_8) = {maxC*W(list(range(8))):.5f}  vs  margin_8={margin8:.5f}  "
      f"=> {'certificate-only CLOSES' if maxC*W(list(range(8)))<margin8 else 'consec near boundary: finite check needed'}")
print("\nDONE.")
