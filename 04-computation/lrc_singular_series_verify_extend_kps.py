#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VERIFY (fresh, independent) the three LRC singular-series claims, then EXTEND via
the near-tight core family.  kind-pasteur-2026-06-15-S3.

Object (THM-501): D(q,S) = #{a in Z/q : v*a mod q not in B_q for all v in S},
B_q = {b : min(b, q-b) <= floor(q/14)} (the dangerous band).  D(q,S)>0 <=> loose
witness at shell q.  L(S) = lim_{q->oo} D(q,S)/q (the singular series; the
covering-density). L>0 ⟹ loose; C'(14) <=> inf L > 0 over primitive mult-of-14 S.

CLAIMS:
 (1) inf L ≈ 0.0237 over primitive mult-of-14 configs; extremizers = {1..12} U {14m}.
 (2) s(t)=sin(pi t/7)/(pi t): sign = sign(sin(pi t/7)), period 14, + on {1..6}, 0 at 7,
     - on {8..13}, 0 at 14; only 7-coprime coeffs contribute.
 (3) adelic: for the evader, D(p^e)/p^e -> L itself for EVERY prime p (not a per-prime
     factor) -> no Euler product (refutes HYP-2503).
EXTENSION: compute L for the near-tight family {1..12}U{14m}, m=1..; estimate the
m->oo limit (the candidate infimum) and whether it is attained or approached.
"""
import sys
from math import gcd, sin, pi
from functools import reduce
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def band_safe(q):
    """return a function/set: a residue r is DANGEROUS if min(r,q-r)<=floor(q/14)."""
    h=q//14
    dang=bytearray(q)
    for r in range(0,h+1): dang[r]=1
    for r in range(q-h,q): dang[r]=1
    dang[0]=1
    return dang

def D(q,S):
    dang=band_safe(q); c=0
    Sr=[v % q for v in S]
    for a in range(q):
        safe=True
        for v in Sr:
            if dang[(v*a)%q]: safe=False; break
        if safe: c+=1
    return c

def L_est(S, Q0=4000, W=400):
    """average D(q,S)/q over q in [Q0, Q0+W]."""
    tot=0.0
    for q in range(Q0,Q0+W):
        tot += D(q,S)/q
    return tot/W

def primitive(S): return reduce(gcd,S)==1

def main():
    # ---- Claim 2: sign structure of s(t) ----
    print("=== CLAIM 2: s(t)=sin(pi t/7)/(pi t) sign structure (period 14) ===")
    def s(t): return 0.0 if t%7==0 else sin(pi*t/7)/(pi*t)
    signs={t: ('+' if s(t)>1e-12 else ('-' if s(t)<-1e-12 else '0')) for t in range(1,15)}
    print("  t:    ", list(range(1,15)))
    print("  sign: ", [signs[t] for t in range(1,15)])
    ok = all(signs[t]=='+' for t in range(1,7)) and signs[7]=='0' and all(signs[t]=='-' for t in range(8,14)) and signs[14]=='0'
    print(f"  + on 1..6, 0 at 7, - on 8..13, 0 at 14:  {ok}   (7-coprime coeffs only contribute)")

    # ---- Claim 1: infimum + extremizers ----
    print("\n=== CLAIM 1: inf L over primitive mult-of-14 configs ===")
    tight = list(range(1,14))  # {1..13}, M=1/14, should give L~0
    print(f"  tight {{1..13}} (M=1/14): L = {L_est(tight):.4f}  (expect ~0)")
    print("  near-tight core family {1..12} U {14m}:")
    fam={}
    for m in range(1,9):
        S=sorted(list(range(1,13))+[14*m])
        if not primitive(S):
            print(f"    m={m}: NOT primitive (skip)"); continue
        Lv=L_est(S); fam[m]=Lv
        print(f"    m={m}: S={{1..12}} U {{{14*m}}}  L = {Lv:.4f}")
    # generic + evader for comparison
    gen=[14,17,19,23,29,31,37,41,43,47,53,59,61]
    evad=sorted([7*k for k in range(1,13)]+[611])
    print(f"  generic coprime: L = {L_est(gen):.4f}  (expect ~ (6/7)^13 = {(6/7)**13:.4f})")
    print(f"  evader 7*{{1..12}}+611: L = {L_est(evad):.4f}  (claim ~0.0293)")

    # ---- Claim 3: adelic test (evader) ----
    print("\n=== CLAIM 3: adelic D(p^e)/p^e -> L for each prime p (evader) ===")
    for p in (2,3,5,7):
        seq=[]
        e=1
        while p**e <= 60000:
            q=p**e; seq.append(D(q,evad)/q); e+=1
        print(f"  p={p}: D(p^e)/p^e for e=1.. = {[f'{x:.4f}' for x in seq]}  -> {seq[-1]:.4f}")
    print(f"  (each prime's sequence should approach L_evader ~0.029, NOT a per-prime factor)")

    # ---- EXTENSION: the m->oo limit of the near-tight family ----
    print("\n=== EXTENSION: trend of L({1..12} U {14m}) as m grows (candidate infimum) ===")
    print("  m, L:", {m:round(v,4) for m,v in fam.items()})
    if len(fam)>=4:
        vals=list(fam.values())
        print(f"  min over computed m: {min(vals):.4f} at m={min(fam,key=fam.get)}")
        print("  (does L decrease toward a limit ~0.0237, or is there an interior min?)")

if __name__=="__main__":
    main()
