#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(n) RECURSION as the number of runners n changes.
kind-pasteur-2026-06-19 (ANGLE: the recursion as n changes).

For LRC(n) (n runners), the singular-series / sector attack uses:
  * lonely floor 1/n  (the conjectured gap).
  * sectors of width 1/s where s = n/2 (even n), the "1/s pivot": a free fast phase
    exists iff the cluster phases leave a gap > 1/s. s = n/2 for n=14 gives s=7.
  * relations among integer co-offsets live mod C = 2n-1.
  * cover set S_s(E) = {x : every sector [j/s,(j+1)/s) hit by some frac(e_i x)},
    e_i in E, 0 in E, |E| = k.  meas(S_s(E)) <= cap_k = min_{|P|=n-1-k} meas(G_P).
  * "consec maximizes" is the crux on the DANGEROUS rows k.

This script computes, FOR A RANGE OF n, the analogous structure:
  - sector count s(n), modulus C=2n-1 and its factorization (= strata count),
  - meas(S_s(consec_k)) and the AP-to-cap margin per cluster size k,
  - which rows are "dangerous" (small margin),
  - tests whether n=14 (C=3^3 prime power) is anomalous.

EXACT via Fraction breakpoints x = a/(s*e).  (s replaces the hardcoded 7.)
"""
import sys
from fractions import Fraction
from functools import reduce
from math import gcd
try:
    from sympy import factorint
    HAVE_SYMPY = True
except Exception:
    HAVE_SYMPY = False

def _factorint(m):
    if HAVE_SYMPY:
        return dict(factorint(m))
    f={}; d=2; x=m
    while d*d<=x:
        while x%d==0:
            f[d]=f.get(d,0)+1; x//=d
        d+=1
    if x>1: f[x]=f.get(x,0)+1
    return f

def N_at(E, x, s):
    """number of sectors among {1..s-1} missed by {frac(e x): e in E} at x (Fraction)."""
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator)
        hit.add((v.numerator*s)//v.denominator)
    return sum(1 for j in range(1,s) if j not in hit)

def dist_p(E, s):
    """exact distribution p_t (Fraction), t=0..s-1, of N over x in [0,1)."""
    E=sorted(set(E))
    bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,s*e+1):
            bps.add(Fraction(a,s*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[Fraction(0)]*s
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        t=N_at(E,(lo+hi)/2,s)
        p[t]+=(hi-lo)
    return p

def meas_S(E, s):
    """meas(S_s(E)) = p_0 = measure of x where all s sectors are hit."""
    return dist_p(E,s)[0]

def consec(k):
    return list(range(k))

def lrc_params(n):
    """sector count s, modulus C, factorization, strata count, speeds n-1."""
    assert n%2==0, "even n only (s=n/2)"
    s=n//2
    C=2*n-1
    fac=_factorint(C)
    return dict(n=n, s=s, speeds=n-1, C=C, fac=fac,
                n_distinct_primes=len(fac),
                max_exp=max(fac.values()),
                is_prime_power=(len(fac)==1),
                is_nontrivial_pp=(len(fac)==1 and max(fac.values())>=2),
                n_strata=sum(fac.values())+1)  # filtration depth for prime power = exp+1

if __name__=="__main__":
    if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
    print("="*78)
    print("LRC(n) RECURSION: sector count, modulus, strata as n changes")
    print("="*78)
    ns=[8,10,12,14,16,18,20]
    print(f"\n{'n':>3} {'s=n/2':>6} {'speeds':>7} {'C=2n-1':>7} {'factor(C)':>14} "
          f"{'#primes':>8} {'maxexp':>7} {'pp?':>5} {'depth':>6}")
    for n in ns:
        P=lrc_params(n)
        facstr="*".join(f"{p}^{e}" if e>1 else f"{p}" for p,e in sorted(P['fac'].items()))
        print(f"{n:>3} {P['s']:>6} {P['speeds']:>7} {P['C']:>7} {facstr:>14} "
              f"{P['n_distinct_primes']:>8} {P['max_exp']:>7} "
              f"{'YES' if P['is_prime_power'] else 'no':>5} {P['n_strata']:>6}")

    print("\n" + "="*78)
    print("Per-n: meas(S_s(consec_k)) for cluster sizes k (the dangerous-row analog)")
    print("Dangerous k for n=14 are 8,9,10 = (speeds-1) down a bit. General: k near n-4..n-1")
    print("="*78)
    for n in ns:
        P=lrc_params(n); s=P['s']; speeds=P['speeds']
        # dangerous cluster sizes: k from about speeds-5 to speeds-1 (mirror k=8,9,10 at speeds=13)
        # at n=14 speeds=13, dangerous k=8,9,10 -> speeds-5..speeds-3 ; also show up to speeds-1.
        klo=max(2, speeds-5); khi=speeds-1
        ks=list(range(klo,khi+1))
        print(f"\n--- n={n}  s={s} sectors  speeds={speeds}  C={P['C']}  k in {ks} ---")
        prev=None
        for k in ks:
            ms=meas_S(consec(k), s)
            # p_{s-1} all-missed analog = 1/(s*(k-1)) (the closed form from THM)
            allm=dist_p(consec(k),s)[s-1]
            print(f"   k={k:2d}: meas(S_s(consec)) = {float(ms):.5f}  ({str(ms)});  "
                  f"p_all-missed={float(allm):.5f}")
