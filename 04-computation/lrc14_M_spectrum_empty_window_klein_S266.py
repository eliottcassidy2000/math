#!/usr/bin/env python3
"""
lrc14_M_spectrum_empty_window_klein_S266.py
===========================================
klein-2026-07-12-S266

VERIFY the consolidated frontier's central claim (opus-S246 HYP-6155 + boxeph-S20 HYP-6150):
the loneliness spectrum M(S) has an EMPTY WINDOW just above the tight value 1/14 — nothing
lands strictly in (1/14, 2/27) — and the next rungs are the Farey fractions r/(13r+1) and
1/13. Equivalently: [M < 2/27 => S is the dilated interval {1..13}] and [DC => M >= 1/13].

M(S) = max_t min_i ||v_i t||, exact via the THM-668 pair-sum ruler (M = p/q, q = v_i+v_j).

This script (all exact rationals):
 (1) confirms the named extremal M-values: interval {1..13}=1/14; the compressed near-dilate
     2*{1..12}u{13}=1/13; boxeph's step-3 AP core =3/37; check they are what canon says.
 (2) enumerates the LOW END of the M-spectrum over primitive 13-subsets of [1..V] for
     V=14..VMAX, records every distinct M <= 0.11, and flags any M in the (1/14, 2/27) window.
 (3) reports the Farey rungs r/(13r+1) and where the observed spectrum sits relative to them.
 (4) the DC (divisor-complete) subset's minimum M in-range (cross-check boxeph's 1/13 floor).
"""
import math
from fractions import Fraction
from itertools import combinations

VMAX = 20          # enumerate primitive 13-subsets of [1..V], V=14..VMAX
LOW  = Fraction(11,100)   # record M <= 0.11

def dist0(r,q): return min(r, q-r)

def exact_M(v, early_stop=None):
    """M(S) via pair-sum rulers.  If early_stop given (Fraction), return as soon as a
       witness proves M >= early_stop (we only care about SMALL M here, so we can bail
       out of a family the moment it's clearly loose)."""
    n=len(v); best=Fraction(0); argq=None
    qs=set()
    for a in range(n):
        for b in range(a,n):
            qs.add(v[a]+v[b])
    for q in sorted(qs):
        mx=0
        for p in range(1,q):
            mn=q
            for vl in v:
                d=dist0((vl*p)%q,q)
                if d<mn:
                    mn=d
                    if mn<=mx: break
            if mn>mx: mx=mn
        val=Fraction(mx,q)
        if val>best:
            best=val; argq=q
            if early_stop is not None and best>=early_stop:
                return best,argq
    return best,argq

def primitive(v):
    g=0
    for x in v: g=math.gcd(g,x)
    return g==1

def divisor_complete(v):
    return all(any(x%d==0 for x in v) for d in range(2,15))

one14=Fraction(1,14); two27=Fraction(2,27); one13=Fraction(1,13)

print("="*72)
print("(1) NAMED EXTREMAL M-VALUES (exact, via pair-sum ruler)")
print("="*72)
named = [
    ("interval {1..13}", list(range(1,14))),
    ("compressed near-dilate 2*{1..12} u {13}", sorted([2*i for i in range(1,13)]+[13])),
    ("boxeph step-3 AP core", [2,3,5,8,9,11,12,13,14,15,17,20,23]),
    ("translated interval {2..14}", list(range(2,15))),
]
for nm,v in named:
    M,q=exact_M(v)
    print(f"  {nm:42s} M={str(M):>8s}={float(M):.5f}  (ruler q={q})  DC={divisor_complete(v)}")
print(f"  reference: 1/14={float(one14):.5f}  2/27={float(two27):.5f}  1/13={float(one13):.5f}")

print()
print("="*72)
print(f"(2) LOW END of the M-spectrum over primitive 13-subsets of [1..V], V=14..{VMAX}")
print("="*72)
spectrum = {}   # M (Fraction) -> (count, example family, anyDC)
window_hits = []
total=0
for V in range(14, VMAX+1):
    # enumerate 13-subsets of [1..V] that CONTAIN V (new max) to avoid recount across V
    rest = range(1, V)
    cnt_V=0
    for combo in combinations(rest, 12):
        v = list(combo)+[V]
        if v[0] != 1:   # primitive fast-screen: interval-like small M needs 1 present;
            # still allow non-1-containing but they are rarely small-M; keep for correctness
            pass
        if not primitive(v):
            continue
        total+=1; cnt_V+=1
        M,_ = exact_M(v, early_stop=LOW+Fraction(1,1000))  # bail once clearly > 0.11
        if M<=LOW:
            key=M
            dc=divisor_complete(v)
            if key not in spectrum:
                spectrum[key]=[0,v,False]
            spectrum[key][0]+=1
            if dc: spectrum[key][2]=True
            if one14 < M < two27:
                window_hits.append((M,v))
    print(f"  V={V:2d}: {cnt_V:>7d} primitive families scanned")

print(f"\n  total primitive families scanned: {total}")
print(f"\n  DISTINCT M-values <= 0.11, sorted (M, count, example, DC?):")
for M in sorted(spectrum):
    c,ex,dc = spectrum[M]
    tag = "  <-- in (1/14,2/27) WINDOW!" if one14<M<two27 else ("  = 1/14 (tight)" if M==one14 else "")
    dctag = " [DC-achievable]" if dc else ""
    print(f"    {str(M):>10s} = {float(M):.5f}  x{c:<6d}  e.g.{ex}{dctag}{tag}")

print()
print("="*72)
print("(3) EMPTY-WINDOW CHECK (1/14, 2/27)")
print("="*72)
if window_hits:
    print(f"  !! {len(window_hits)} families with M strictly in (1/14, 2/27) -- WINDOW NOT EMPTY:")
    for M,v in window_hits[:10]:
        print(f"     M={M}={float(M):.5f}  {v}")
else:
    print(f"  EMPTY: no primitive 13-set in [1..{VMAX}] has M strictly in (1/14, 2/27).")
    print(f"  => the spectrum jumps from 1/14 directly to the next rung >= 2/27 (verified in-range).")

print()
print("="*72)
print("(4) FAREY RUNGS r/(13r+1) vs observed, and the DC floor in-range")
print("="*72)
print("  Farey rungs r/(13r+1):", [f"{r}/{13*r+1}={float(Fraction(r,13*r+1)):.5f}" for r in range(1,5)])
dc_ms = [M for M in spectrum if spectrum[M][2]]
if dc_ms:
    print(f"  min DC M in-range = {min(dc_ms)} = {float(min(dc_ms)):.5f}  (boxeph floor: 1/13={float(one13):.5f})")
else:
    print("  no DC families with M<=0.11 in-range (DC needs Vmax>=14; compressed extremal at Vmax=26 out of range)")
print("\ndone.")
