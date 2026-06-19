#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_mu17lean_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL, LEAN)
Fast targeted hunt for mu_1/7(E) < mu_1/7(consec_k), k=8..13.
Exhaustive up to a MODEST spread (where C(s-1,k-2) feasible) + heavy structured/large-spread.
Float prefilter (fast), exact Fraction confirmation only on near-ties.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(424242)
ONE7 = F(1,7); ONE7f = 1.0/7.0

# fast float maxgap-measure via collision breakpoints + midpoint sampling at high resolution.
# For a RIGOROUS check we use exact engine; for the HUNT we use a fast exact engine but only on
# candidates. To keep speed, compute mu_1/7 exactly with the cell engine (it's fine for k<=13 when
# spread is modest). We add a float prescreen for the random shower.

def mu17_exact(E):
    E = sorted(set(E)); k=len(E)
    if k==1: return F(1)
    bp={F(0),F(1)}
    for i in range(k):
        for j in range(i+1,k):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bps=sorted(b for b in bp if F(0)<=b<=F(1))
    tot=F(0)
    for a,b in zip(bps,bps[1:]):
        if a==b: continue
        mid=(a+b)/2
        pts=[]
        for e in E:
            val=e*mid; n=val-(val%1); pts.append((val-n,e,n))
        pts.sort(key=lambda t:t[0]); m=len(pts); ivs=[]
        for i in range(m):
            (_,ei,ni)=pts[i]
            if i<m-1: (_,ej,nj)=pts[i+1]; al=ej-ei; be=ni-nj
            else: (_,e0,n0)=pts[0]; al=e0-ei; be=ni-n0+1
            al=F(al); be=F(be)
            if al==0:
                if be>ONE7: ivs.append((a,b))
            else:
                xs=(ONE7-be)/al
                if al>0: lo,hi=max(a,xs),b
                else: lo,hi=a,min(b,xs)
                if lo<hi: ivs.append((lo,hi))
        if not ivs: continue
        ivs.sort(); clo,chi=ivs[0]
        for lo,hi in ivs[1:]:
            if lo<=chi: chi=max(chi,hi)
            else: tot+=chi-clo; clo,chi=lo,hi
        tot+=chi-clo
    return tot

def mu17_float(E, N=20000):
    E=sorted(set(E)); cnt=0
    for t in range(N):
        x=(t+0.5)/N
        pts=sorted((e*x)%1.0 for e in E)
        mg=0.0
        for i in range(len(pts)-1):
            g=pts[i+1]-pts[i]
            if g>mg: mg=g
        w=pts[0]+1.0-pts[-1]
        if w>mg: mg=w
        if mg>ONE7f: cnt+=1
    return cnt/N

def gcd1(E): return reduce(gcd,E)==1

consec={k:mu17_exact(list(range(k))) for k in range(8,14)}
print("consec mu_1/7:", {k:str(consec[k]) for k in consec})

print("\n=== EXHAUSTIVE (modest spread, exact) ===")
CAP=12000
exh_ok=True
for k in range(8,14):
    base=consec[k]; best=base; bestE=None
    s=k-1; smax=k-1
    while comb(s-1,k-2)<=CAP: smax=s; s+=1
    cnt=0
    for s in range(k-1,smax+1):
        for interior in itertools.combinations(range(1,s),k-2):
            E=(0,)+interior+(s,)
            if not gcd1(E): continue
            cnt+=1
            m=mu17_exact(list(E))
            if m<best: best=m; bestE=E
    if best<base: exh_ok=False
    print(f"  k={k}: spread<= {smax} ({cnt} sets) min={best}={float(best):.5f} consec={float(base):.5f} "
          f"{'OK' if best>=base else 'VIOLATION '+str(bestE)}")
print(f"  exhaustive-modest: consec-is-min = {exh_ok}")

print("\n=== STRUCTURED + LARGE-SPREAD (float prescreen, exact confirm on ties) ===")
struct_ok=True
for k in range(8,14):
    basef=float(consec[k]); base=consec[k]
    cands=set()
    # perforated APs over {0..k-1+extra}
    for extra in range(1,10):
        full=list(range(k+extra))
        for drops in itertools.combinations(range(1,k+extra),extra):
            E=tuple(x for x in full if x not in drops)
            if len(E)==k and E[0]==0: cands.add(E)
    # common-factor / subtorus mixes
    for c in (2,3,4,5,6,7):
        for tlen in range(2,k):
            A=[c*i for i in range(tlen)]
            rem=k-tlen
            B=list(range(1,rem+1))
            E=tuple(sorted(set([0]+A[1:]+B)))
            if len(E)==k and 0 in E: cands.add(E)
    # two-scale
    for t in range(2,k):
        rem=k-(t+1)
        if rem<1: continue
        for M in (k,k+2,2*k,3*k,5*k,10*k,30*k,100*k):
            E=tuple(sorted(set(list(range(t+1))+list(range(M,M+rem)))))
            if len(E)==k and 0 in E: cands.add(E)
    # big random shower
    for _ in range(15000):
        sp=random.choice([k-1,k,k+1,k+2,k+3,k+5,k+8,2*k,3*k,5*k,10*k,20*k])
        if sp<k-1: continue
        body=random.sample(range(1,sp+1),min(k-1,sp))
        E=tuple(sorted(set([0]+body)))
        if len(E)==k: cands.add(E)
    viol=[]
    for E in cands:
        if len(set(E))!=k or 0 not in E: continue
        mf=mu17_float(list(E),12000)
        if mf < basef + 0.01:   # near or below -> exact confirm
            me=mu17_exact(list(E))
            if me<base: viol.append((me,E))
    if viol:
        struct_ok=False; viol.sort()
        print(f"  k={k}: *** {len(viol)} VIOLATIONS *** lowest mu={viol[0][0]}={float(viol[0][0]):.5f} at {viol[0][1]}")
    else:
        print(f"  k={k}: {len(cands)} cands, no E beats consec {basef:.5f}")
print(f"  structured/large: consec-is-min = {struct_ok}")
print(f"\nVERDICT 1/7 consecutive-minimizes (lean): exhaustive={exh_ok} structured={struct_ok} "
      f"=> {exh_ok and struct_ok}")
