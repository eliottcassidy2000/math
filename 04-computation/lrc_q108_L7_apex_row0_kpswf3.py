#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_L7_apex_row0_kpswf3.py  (kind-pasteur 2026-06-21, THREAD A proof of step D)

Closed form of ROW 0 of mu_{p,q}, and the proof that row-0 is flat <=> 7|p (7 nmid q).

Row 0 = { v : sector(qv)=0 } = { v in union_{k=0}^{q-1} [k/q, k/q + 1/(7q)) } (q bands,
each length 1/(7q), one per period of qv).  We must distribute, across these q bands,
the p-sector j = floor(7 p v) mod 7.

On band k: v = k/q + w/(7q), w in [0,1).  Then
   7 p v = 7 p k / q + p w / q  =  7pk/q + (p/q) w.
   j = floor(7 p v) mod 7 = floor( 7pk/q + (p/q) w ) mod 7.
As w sweeps [0,1), (p/q)w sweeps [0, p/q) -- length p/q, so it crosses floor-boundaries
p times (gives p unit steps) ... the band's mass 1/(7q) is split among the j-values it visits.

CLAIM proven here numerically + the arithmetic reason:
  row0(j) = (1/(7q)) * Leb{ w in [0,1) : floor( 7pk/q + (p/q) w ) ≡ j  (mod 7) } summed over k.
The j visited by band k starts at  J_k := floor(7pk/q) mod 7  and the offset of the start
within its sector is determined by frac(7pk/q). Summing the q bands, row0 is flat
<=> the q starting offsets { 7pk/q mod 7 : k=0..q-1 } equidistribute mod 7 with the right
weights, which (since gcd(p,q)=1) happens EXACTLY when 7 | p (the slope is a multiple of a
full sector per band) -- otherwise the q bands land on a proper coset pattern mod 7.

We verify the EXACT closed form and the divisibility equivalence.
"""
from fractions import Fraction as Fr
from math import gcd
P = 7
def sector(yf): return int(P * yf)

def mu_full(p, q):
    bp = {Fr(0), Fr(1)}
    for f in (p, q):
        for t in range(0, P * f): bp.add(Fr(t, P * f))
    vs = sorted(bp); cell = {}
    for a, b in zip(vs, vs[1:]):
        mid = (a + b) / 2
        cell[(sector((q*mid)%1), sector((p*mid)%1))] = \
            cell.get((sector((q*mid)%1), sector((p*mid)%1)), Fr(0)) + (b - a)
    return cell

def row0_direct(p, q):
    """row0(j) by the band decomposition (independent recomputation)."""
    out = {j: Fr(0) for j in range(P)}
    for k in range(q):
        # band: v in [k/q, k/q + 1/(7q)); param v = k/q + w/(7q), w in [0,1)
        # breakpoints of j = floor(7 p v) inside the band:
        lo = Fr(k, q); hi = Fr(k, q) + Fr(1, 7*q)
        # collect v-breakpoints where 7pv crosses an integer
        bps = {lo, hi}
        # 7 p v integer => v = t/(7p)
        t0 = (7*p*lo).__ceil__()
        t1 = (7*p*hi).__floor__()
        for t in range(t0, t1+1):
            vv = Fr(t, 7*p)
            if lo <= vv <= hi: bps.add(vv)
        bl = sorted(bps)
        for a, b in zip(bl, bl[1:]):
            mid = (a+b)/2
            j = sector((p*mid) % 1)
            out[j] += (b - a)
    return [out[j] for j in range(P)]

# verify row0_direct matches mu_full's row 0, and the flat<=>7|p law
ok_match = True; ok_law = True
for q in range(1, 40):
    if q % P == 0: continue
    for p in range(1, 40):
        if gcd(p, q) != 1: continue
        cell = mu_full(p, q)
        r0a = [cell.get((0, j), Fr(0)) for j in range(P)]
        r0b = row0_direct(p, q)
        if r0a != r0b: ok_match = False
        flat = all(x == Fr(1, 49) for x in r0b)
        if flat != (p % P == 0): ok_law = False
print("row0 band-decomposition matches mu_full row 0:", "YES" if ok_match else "NO")
print("row0 flat <=> 7|p (given 7 nmid q):", "HOLDS" if ok_law else "FAILS")

# Show the explicit row0 vectors for a few cases (times 49 for readability)
print("\n49 * row0(j) for sample (p,q):")
for (p, q) in [(7,3),(2,3),(3,5),(1,4),(14,5),(5,3),(8,3)]:
    if gcd(p,q)!=1 or q%7==0: continue
    r0 = row0_direct(p, q)
    print(f"  p={p:>2} q={q}: {[str(x*49) for x in r0]}  (7|p={p%7==0}; flat={all(x==Fr(1,49) for x in r0)})")
