#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(n) RECURSION part 2: the AP-to-cap MARGIN across n.
kind-pasteur-2026-06-19.

cap_k(n) = min_{|P|=n-1-k} meas(G_P),  G_P = {x in [0,1): ||p x|| >= 1/n  for all p in P},
  ||y|| = distance to nearest integer.  (THM-530: for n=14 the mins over |P|=5 give 2243/5880.)
The crux margin for the DANGEROUS row k is
  margin_k(n) = cap_k(n) - meas(S_s(consec_k))      (must be > 0; consec extremal closes LRC(n)).

meas(G_P): G_P is a union of arcs; exact via breakpoints x = (a +/- 1/n)/p.
We minimize meas(G_P) over P (small integer speed sets with the right size) to get cap_k.
Since the true min over ALL valid P is expensive, we min over a representative family:
consecutive-type and small-spread P (THM-530: min is attained at structured small P).

Outputs per n: for each dangerous k, meas(S_s(consec_k)), best (smallest) cap_k found,
the margin, and whether the row is binding (small margin).  Then correlate margin with
the factorization of C=2n-1.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
try:
    from sympy import factorint
    HAVE_SYMPY=True
except Exception:
    HAVE_SYMPY=False

def _factorint(m):
    if HAVE_SYMPY: return dict(factorint(m))
    f={}; d=2; x=m
    while d*d<=x:
        while x%d==0: f[d]=f.get(d,0)+1; x//=d
        d+=1
    if x>1: f[x]=f.get(x,0)+1
    return f

# ---- meas(S_s(E)): all-sectors-hit measure (the S7 analog) ----
def meas_S(E, s):
    E=sorted(set(E))
    bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,s*e+1):
            bps.add(Fraction(a,s*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p0=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        x=(lo+hi)/2
        hit=set()
        for e in E:
            v=e*x; v=v-(v.numerator//v.denominator)
            hit.add((v.numerator*s)//v.denominator)
        if all(j in hit for j in range(s)):
            p0+=(hi-lo)
    return p0

# ---- meas(G_P) = meas{x: ||p x|| >= 1/n for all p in P} ----
def meas_GP(P, n):
    P=sorted(set(p for p in P if p>0))
    if not P: return Fraction(1)
    # breakpoints where ||p x|| crosses 1/n: p x = a +/- 1/n  -> x=(a*n +/- 1)/(p*n)
    bps={Fraction(0),Fraction(1)}
    for p in P:
        for a in range(0,p+1):
            for off in (-1,1):
                x=Fraction(a*n+off, p*n)
                if 0<=x<=1: bps.add(x)
    bps=sorted(bps)
    tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        x=(lo+hi)/2
        ok=True
        for p in P:
            v=p*x; v=v-(v.numerator//v.denominator)  # frac(p x) in [0,1)
            d=min(v,1-v)  # ||p x||
            if d<Fraction(1,n): ok=False; break
        if ok: tot+=(hi-lo)
    return tot

def best_cap(n, psz, max_speed=None):
    """min meas(G_P) over structured P with |P|=psz. Search consecutive-ish small sets."""
    if max_speed is None:
        max_speed=psz+5
    best=None; bestP=None
    # candidate speeds 1..max_speed, choose psz of them; prune by count
    cand=list(range(1,max_speed+1))
    # full combos can blow up; cap the search
    from math import comb
    if comb(len(cand),psz)>20000:
        # restrict to windows
        cand=list(range(1,psz+4))
    for P in itertools.combinations(cand,psz):
        m=meas_GP(P,n)
        if best is None or m<best:
            best=m; bestP=P
    return best,bestP

if __name__=="__main__":
    if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
    print("="*92)
    print("LRC(n) AP-to-cap MARGIN across n   (cap_k = min_|P|=n-1-k meas(G_P); margin=cap-meas_S(consec))")
    print("="*92)
    ns=[8,10,12,14,16,18,20]
    rows=[]
    for n in ns:
        s=n//2; speeds=n-1; C=2*n-1; fac=_factorint(C)
        facstr="*".join(f"{p}^{e}" if e>1 else f"{p}" for p,e in sorted(fac.items()))
        maxexp=max(fac.values()); npr=len(fac)
        print(f"\n--- n={n}  s={s}  speeds={speeds}  C={C}={facstr}  (#primes={npr}, maxexp={maxexp}) ---")
        # dangerous rows: mirror n=14 (k=8,9,10 -> 13-k = |P| = 5,4,3). Use |P|=3,4,5 => k=speeds-3..speeds-5
        worst_margin=None; worst_k=None
        for psz in (5,4,3):
            k=speeds-psz
            if k<2: continue
            ms=meas_S(list(range(k)), s)
            cap,Pbest=best_cap(n,psz)
            margin=cap-ms
            tag=" <<< BINDING" if margin< Fraction(2,100) else ""
            print(f"   k={k:2d} (|P|={psz}): meas_S(consec)={float(ms):.5f}  cap={float(cap):.5f} ({str(cap)} P={Pbest})  margin={float(margin):+.5f}{tag}")
            if worst_margin is None or margin<worst_margin:
                worst_margin=margin; worst_k=k
        rows.append((n,C,facstr,npr,maxexp,worst_k,float(worst_margin)))

    print("\n" + "="*92)
    print("SUMMARY: tightest (binding) margin per n vs factorization of C=2n-1")
    print("="*92)
    print(f"{'n':>3} {'C':>4} {'factor':>10} {'#pr':>4} {'maxexp':>7} {'tight k':>8} {'min margin':>11}")
    for (n,C,facstr,npr,maxexp,wk,wm) in rows:
        print(f"{n:>3} {C:>4} {facstr:>10} {npr:>4} {maxexp:>7} {wk:>8} {wm:>+11.5f}")
