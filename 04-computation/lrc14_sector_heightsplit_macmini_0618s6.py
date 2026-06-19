#!/usr/bin/env python3
"""
lrc14_sector_heightsplit_macmini_0618s6.py  (mac-mini-2026-06-18-S6)

The RELATION-HEIGHT SPLIT for the seven-sector cover meas(S7(E)) <= cap_k (HYP-2603 => HYP-2602).
meas(S7(E)) = M7(k) + corr(E),  M7(k)=Sum_{T subset {1..6}}(-1)^|T|(1-|T|/7)^{k-1} (the iid 1-dim main),
corr(E) = the relation-lattice correction (the singular-series tail, like HYP-2600a/2601).

GOAL: show (1) M7(k) << cap_k so the certificate margin_k=cap_k-M7(k) is large;
      (2) corr(E) is SMALL for high-relation-height (dissociated) E -> certified;
      (3) corr(E) is maximized by CONSECUTIVE (the binding low-height case), still <= margin_k;
      (4) the low-height residual (uncertified) is the AP-rich family, bounded modulo scale.
Output: M7(k), cap_k, margin_k (k=8..13); for shapes: meas(S7), corr, min-triple-height H_min,
        whether consec maximizes meas(S7).
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

def M7(k):
    return sum(F((-1)**t*comb(6,t))*F(7-t,7)**(k-1) for t in range(7))

# m_k = cap_k = min over |P|=13-k of meas(G_P)
H=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H/u)%1; b=(c+H/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mg(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    if not P: return F(1)
    dz=mg([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
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

def min_triple_height(E):
    """min |a*b*c| over primitive 3-term relations a*e_i+b*e_j+c*e_l=0 among offsets."""
    E=sorted(set(E)); best=None
    for i,j,l in itertools.combinations(range(len(E)),3):
        a,b,c = E[j]-E[l], E[l]-E[i], E[i]-E[j]   # a*e_i+b*e_j+c*e_l=0 identically
        g=gcd(gcd(abs(a),abs(b)),abs(c)) or 1
        a,b,c=a//g,b//g,c//g
        h=abs(a*b*c)
        if h>0 and (best is None or h<best): best=h
    return best

print("="*78)
print("M7(k), cap_k, margin_k = cap_k - M7(k)  (the certificate budget)")
print("="*78)
for k in range(8,14):
    m=M7(k); c=cap(k);
    print(f"  k={k:2d}: M7={float(m):.5f}  cap_k={float(c):.5f}  margin_k={float(c-m):.5f}")

print("\n"+"="*78)
print("per-shape (k=8): meas(S7), corr=meas(S7)-M7, min-triple-height, meas(S7)<=cap_8?")
print("="*78)
M78=M7(8); cap8=cap(8)
shapes=[
  ("consec {0..7}", list(range(8))),
  ("perf {0,2..7,9}", [0,2,3,4,5,6,7,9]),
  ("3-dilate consec", [3*i for i in range(8)]),
  ("dissoc 1,3,7,15..", [0,1,3,7,15,31,63,127]),
  ("Sidon", [0,1,3,7,12,20,30,44]),
  ("emb-AP {0..3,40..43}", [0,1,2,3,40,41,42,43]),
  ("generic", [0,5,13,27,41,58,79,97]),
]
print(f"  M7(8)={float(M78):.5f}  cap_8={float(cap8):.5f}  margin_8={float(cap8-M78):.5f}")
for name,E in shapes:
    s7=measS7(E); corr=s7-M78; Hm=min_triple_height(E)
    print(f"  {name:<22} meas(S7)={float(s7):.5f}  corr={float(corr):+.5f}  Hmin={Hm}  "
          f"{'OK' if s7<=cap8 else 'VIOL'}")

print("\n"+"="*78)
print("Does CONSECUTIVE maximize meas(S7)? exhaustive k=8 over shapes spread<=12")
print("="*78)
cons=measS7(list(range(8))); mx=(cons,'consec'); viol=0; checked=0
for rest in itertools.combinations(range(1,13),7):
    E=[0]+list(rest); checked+=1; s=measS7(E)
    if s>mx[0]: mx=(s,E)
    if s>cons: viol+=1
print(f"  meas(S7,consec)={float(cons):.5f}; checked {checked}; shapes with S7>consec: {viol}; "
      f"max={float(mx[0]):.5f} {'at consec' if mx[1]=='consec' else 'at '+str(mx[1])}")
print("\nDONE.")
