#!/usr/bin/env python3
"""
lrc14_sector_certificate_macmini_0618s6.py  (mac-mini-2026-06-18-S6)

Build the high-height certificate for meas(S7(E)) <= cap_k and pin the low-height residual.
 (1) finite check: meas(S7(consec_k)) <= cap_k for k=8..13 (is the extremal within budget?).
 (2) consec maximizes meas(S7): exhaustive k=9,10 (spread<=k+3).
 (3) a CONCRETE pairwise correction bound: corr(E) is driven by PAIRS (e_a,e_b) of offsets
     whose orbit overlaps in the sectors beyond independence. Compute the exact pairwise
     excess  P2(E)=sum_{a<b} [meas{x: e_a x, e_b x in same/related sectors} structure].
     Proxy: short-relation weight W(E)=sum over primitive triples of 1/height. Correlate
     with corr(E); the certificate is "few short relations => corr small => certified".
 (4) low-height residual: which shapes have meas(S7) within 0.9*cap_k (near binding)?
     characterize them (AP-rich, bounded mod scale).
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

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

def short_rel_weight(E, Hmax=200):
    """sum over primitive 3-term relations (height<=Hmax) of 1/height — a low-height richness proxy."""
    E=sorted(set(E)); w=F(0)
    for i,j,l in itertools.combinations(range(len(E)),3):
        a,b,c=E[j]-E[l],E[l]-E[i],E[i]-E[j]
        g=gcd(gcd(abs(a),abs(b)),abs(c)) or 1
        h=abs((a//g)*(b//g)*(c//g))
        if 0<h<=Hmax: w+=F(1,h)
    return w

print("="*74)
print("(1) finite check: meas(S7(consec_k)) <= cap_k ?  (the extremal within budget)")
print("="*74)
for k in range(8,14):
    s=measS7(list(range(k))); c=cap(k)
    print(f"  k={k:2d}: meas(S7,consec)={float(s):.5f}  cap_k={float(c):.5f}  "
          f"slack={float(c-s):+.5f}  {'OK' if s<=c else '*** VIOLATES ***'}")

print("\n"+"="*74)
print("(2) consec maximizes meas(S7): exhaustive k=9,10")
print("="*74)
for k,sp in [(9,12),(10,13)]:
    cons=measS7(list(range(k))); mx=(cons,'consec'); viol=0; checked=0
    for rest in itertools.combinations(range(1,sp+1),k-1):
        E=[0]+list(rest); checked+=1; s=measS7(E)
        if s>mx[0]: mx=(s,E)
        if s>cons: viol+=1
    print(f"  k={k}: meas(S7,consec)={float(cons):.5f}; checked {checked}; S7>consec: {viol}; "
          f"max {'at consec' if mx[1]=='consec' else float(mx[0])}")

print("\n"+"="*74)
print("(3) corr(E) vs short-relation weight W(E) (k=8): certificate = W small => corr small")
print("="*74)
M78=M7(8)
shapes=[("consec",list(range(8))),("perf",[0,2,3,4,5,6,7,9]),
        ("dissoc",[0,1,3,7,15,31,63,127]),("Sidon",[0,1,3,7,12,20,30,44]),
        ("emb-AP",[0,1,2,3,40,41,42,43]),("generic",[0,5,13,27,41,58,79,97]),
        ("2blocks",[0,1,2,20,21,22,40,41])]
for name,E in shapes:
    corr=measS7(E)-M78; W=short_rel_weight(E)
    print(f"  {name:<10} corr={float(corr):+.5f}  W(E)={float(W):.4f}  (#offsets pairs etc.)")

print("\n"+"="*74)
print("(4) low-height residual: shapes with meas(S7) >= 0.85*cap_8 (near-binding), k=8")
print("="*74)
cap8=cap(8); thresh=cap8*F(85,100); near=[]
for rest in itertools.combinations(range(1,13),7):
    E=[0]+list(rest); s=measS7(E)
    if s>=thresh: near.append((float(s),E))
near.sort(reverse=True)
print(f"  threshold 0.85*cap_8={float(thresh):.4f}; #near-binding shapes: {len(near)}")
for s,E in near[:8]: print(f"     meas(S7)={s:.4f}  E={E}")
print("  (these AP-rich shapes are the low-height residual = finite mod scale/translation)")
print("\nDONE.")
