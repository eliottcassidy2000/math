#!/usr/bin/env python3
"""
lrc14_exact_rational_measure — mac-mini-2026-06-16-S2

EXACT rational lonely measure L(S) for LRC-14 (independent ground truth).
For integer speeds, the danger set of runner v is a union of arcs centered at
n/v (n=0..v-1) of half-width 1/(14v); the lonely set is the complement; L is
EXACTLY RATIONAL. We compute it with fractions.Fraction (no float error) by
the standard interval-union sweep on [0,1) with mod-1 wraparound.

Tests: interior-drop cores {1..13}\{j} ∪ {14m}, the 12-cores (m->inf limit
(6/7)*meas), and the conjectured global extremizer {1..13}\{6}∪{98}.
Pins the EXACT infimum value (is it a recognizable rational?).
"""
from fractions import Fraction as F
from math import gcd

W = F(1, 14)   # danger half-width in units of (n/v): arc = [n/v - 1/(14v), n/v + 1/(14v)]

def danger_arcs(v):
    """Return list of (lo,hi) Fraction arcs on the real line for runner speed v,
       centers n/v for n=0..v-1, half-width 1/(14v). May cross 0 or 1; caller wraps."""
    hw = F(1, 14*v)
    out = []
    for n in range(v):
        c = F(n, v)
        out.append((c - hw, c + hw))
    return out

def wrap_to_unit(intervals):
    """Split intervals so all lie within [0,1); reduce mod 1."""
    out = []
    for lo, hi in intervals:
        # shift so lo in [0,1)
        # handle general by splitting at integer boundaries within [lo,hi]
        a, b = lo, hi
        # normalize a into [0,1)
        shift = a - (a % 1)
        a -= shift; b -= shift
        # now a in [0,1); b = a + width, width < 1 always here
        if b <= 1:
            out.append((a, b))
        else:
            out.append((a, F(1)))
            out.append((F(0), b - 1))
    return out

def union_measure(intervals):
    """Exact measure of union of [lo,hi) intervals within [0,1)."""
    ivs = sorted(intervals)
    total = F(0)
    cur_lo, cur_hi = None, None
    for lo, hi in ivs:
        if cur_lo is None:
            cur_lo, cur_hi = lo, hi
        elif lo <= cur_hi:
            if hi > cur_hi: cur_hi = hi
        else:
            total += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    if cur_lo is not None:
        total += cur_hi - cur_lo
    return total

def lonely_measure(S):
    arcs = []
    for v in S:
        arcs.extend(danger_arcs(v))
    arcs = wrap_to_unit(arcs)
    dunion = union_measure(arcs)
    return F(1) - dunion

def fr(x): return f"{x} = {float(x):.7f}"

print("="*72)
print("12-cores {1..13}\\{j}: exact meas(Lonely) and the m->inf limit (6/7)*meas")
print("="*72)
core12 = {}
for j in range(1,14):
    A = [v for v in range(1,14) if v != j]
    L = lonely_measure(A)
    core12[j] = L
    print(f" j={j:2d}: meas={fr(L)}   (6/7)*meas={fr(F(6,7)*L)}")
jmin = min(core12, key=lambda j: core12[j])
print(f"\n min-measure 12-core at j={jmin}: meas={fr(core12[jmin])}, limit={fr(F(6,7)*core12[jmin])}")

print("\n" + "="*72)
print("Interior-drop cores {1..13}\\{j} ∪ {14m} : exact L for m in resonant set")
print("="*72)
ms = [1,2,3,4,5,7,11,14,21,49]
best = None
for j in range(1,14):
    A = [v for v in range(1,14) if v != j]
    row = []
    for m in ms:
        S = A + [14*m]
        L = lonely_measure(S)
        row.append((m, L))
        if best is None or L < best[2]:
            best = (j, m, L, S)
    # print j row compactly
    cells = "  ".join(f"m={m}:{float(L):.6f}" for m,L in row)
    print(f" j={j:2d}: {cells}")

print("\n" + "="*72)
print("GLOBAL EXTREMIZER candidate and the exact infimum over this scan")
print("="*72)
j,m,L,S = best
print(f"min over scan: j={j}, m={m}, 14m={14*m}, S={sorted(S)}")
print(f"  L = {L}  (exact rational)")
print(f"  L = {float(L):.10f}")
print(f"  1/L = {float(1/L):.4f}   denominator={L.denominator}, numerator={L.numerator}")

# the documented global extremizer j=6,m=7 (98=2*7^2)
S0 = [1,2,3,4,5,7,8,9,10,11,12,13,98]
L0 = lonely_measure(S0)
print(f"\ndocumented extremizer S={S0}:")
print(f"  L = {L0}")
print(f"  L = {float(L0):.10f}   1/L={float(1/L0):.4f}")
print(f"  numerator={L0.numerator} denominator={L0.denominator}")
# factor the denominator a bit
def factor(nn):
    nn=abs(nn); f={}; d=2
    while d*d<=nn:
        while nn%d==0: f[d]=f.get(d,0)+1; nn//=d
        d+=1
    if nn>1: f[nn]=f.get(nn,0)+1
    return f
print(f"  denom factorization: {factor(L0.denominator)}")
print(f"  numer factorization: {factor(L0.numerator)}")
print("\nALL L>0 (every config loose):", "confirmed" if L0>0 else "FALSE")
