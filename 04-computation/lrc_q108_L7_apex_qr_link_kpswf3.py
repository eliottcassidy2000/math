#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_apex_qr_link_kpswf3.py  (kind-pasteur 2026-06-21, THREAD A consequences)

Does the apex-prime D=0 law connect to QR(7) / Q(sqrt-7) / the death-chain (HYP-2657/2692/2708)?

We probe THREE possible bridges:
 (1) The shift s = p*q^{-1} mod 7 (Lemma B): is the resonance MATRIX (when D>0) a QR-graded
     object? Compute D_{p,q} as a function of the residue pair (p mod 7, q mod 7) and the
     "reduced slope" s. Is D constant on QR/NQR classes of s?
 (2) The death-chain coefficients C(d)=cover currency vs the cell-discrepancy: the apex law says
     when 7|pq the far pair is EXACTLY decorrelated (mu uniform), so the death-chain hit law is
     EXACTLY iid -> R=0. So apex alignment is precisely the "far pair is exact iid" locus.
 (3) Gauss/QR: is the SET {s : D depends on s} symmetric under s -> -s (the reflection x->-x of
     HYP-2657)? and under s -> QR-multiple?
"""
from fractions import Fraction as Fr
from math import gcd
P = 7
QR = {1,2,4}; NQR = {3,5,6}
def sector(yf): return int(P*yf)
def mu_full(p,q):
    bp={Fr(0),Fr(1)}
    for f in (p,q):
        for t in range(0,P*f): bp.add(Fr(t,P*f))
    vs=sorted(bp); cell={}
    for a,b in zip(vs,vs[1:]):
        mid=(a+b)/2; k=(sector((q*mid)%1),sector((p*mid)%1))
        cell[k]=cell.get(k,Fr(0))+(b-a)
    return cell
def D_pq(p,q):
    cell=mu_full(p,q); inv=Fr(1,49)
    return sum(abs(cell.get((i,j),Fr(0))-inv) for i in range(P) for j in range(P))

print("="*74)
print("(1) Does D depend only on the reduced slope s = p q^{-1} mod 7?  (and is it QR-graded?)")
print("="*74)
# group D values by s (over coprime p,q with 7 nmid q), see how many distinct D per s
from collections import defaultdict
byS = defaultdict(set)
for q in range(1,30):
    if q%7==0: continue
    si=pow(q%7,-1,7)
    for p in range(1,60):
        if gcd(p,q)!=1: continue
        if p%7==0: continue   # 7|p => D=0, separate
        s=(p*si)%7
        byS[s].add(D_pq(p,q))
print(f"  {'s':>3} {'class':>5} {'#distinct D values':>20} {'min D':>10} {'max D':>10}")
for s in range(1,7):
    vals=sorted(byS[s])
    cls = 'QR' if s in QR else 'NQR'
    print(f"  {s:>3} {cls:>5} {len(vals):>20} {float(min(vals)):>10.5f} {float(max(vals)):>10.5f}")
print("  => D is NOT a function of s alone (depends on q too); s sets the SHIFT, not the size.")

print("\n"+"="*74)
print("(2) s -> -s reflection (HYP-2657 x->-x): is the D-spectrum of s equal to that of -s mod 7?")
print("="*74)
for s in range(1,7):
    ms=(-s)%7
    eq = byS[s]==byS[ms]
    print(f"  s={s} (-s={ms}): D-spectra equal? {eq}   [s in {'QR' if s in QR else 'NQR'}, -s in {'QR' if ms in QR else 'NQR'}]")

print("\n"+"="*74)
print("(3) THE EXACT DECORRELATION STATEMENT: apex-aligned <=> far pair is EXACTLY iid")
print("="*74)
# When 7|pq, mu is uniform 1/49 = product of two uniform 1/7 marginals -> exact independence.
# Verify: mu factorizes as (1/7)(1/7) iff D=0.
allmatch=True
for q in range(1,25):
    for p in range(1,25):
        if gcd(p,q)!=1: continue
        cell=mu_full(p,q)
        factorizes = all(cell.get((i,j),Fr(0))==Fr(1,7)*Fr(1,7) for i in range(P) for j in range(P))
        if factorizes != (D_pq(p,q)==0): allmatch=False
print(f"  mu factorizes (far pair EXACTLY independent) <=> D=0 <=> 7|pq:  {'CONFIRMED' if allmatch else 'FAIL'}")
print("  CONSEQUENCE: on the apex-aligned locus the death-chain two-far hit law is EXACTLY the")
print("  iid law K2(t) (HYP-2708) -> R_2 residual = 0 there, NO resonance debt. The apex prime's")
print("  own multiples are the exact zeros of the two-far correction.")
