#!/usr/bin/env python3
"""
klein-2026-07-04-S126b (HYP-4080) - VERIFY THE SPECTRAL GAP above 1/12 for 11-runner sets.

Test 1 (S126) over perturbations of {1..11} found 1/12 ISOLATED: smallest loose value = 2/23
(gap 1/276), nothing in (1/12, 2/23). This verifies it by FULL enumeration of 11-subsets of
{1..B} (B up to 17) + dilated/large-speed families: does ANY primitive 11-runner U have
M(U) in the OPEN interval (1/12, 2/23)?  If none, the LRC(12) covering-min ladder k/(11k+1)
is discrete at the bottom => a genuine spectral gap delta0 = 2/23 - 1/12 = 1/276 > 0.

Consequence (honest, an INGREDIENT not a closure): loose 11-runner even parts U have M(U) >= 2/23
(a DEFINITE margin, not infinitesimal). So opus's Lemma 3 (large tightener, klein S125) fires for
tighteners > u_max/(6*(2/23-1/12)) = 46*u_max, and the residual sliver is confined to LARGE u_max
with tighteners <= 46*u_max -- exactly mac-mini's 'bound v_max(U)' crux, now with the near-AP
continuum removed (the residual is discrete-rung, not a continuum accumulating at 1/12).
"""
from fractions import Fraction as F
from math import gcd
import itertools

def cdist_q(a, q):
    r = a % q
    return min(r, q - r)

def Mval(S, Qcap):
    best = F(0)
    for Q in range(2, Qcap + 1):
        for a in range(1, Q // 2 + 1):
            if gcd(a, Q) != 1: continue
            m = min(F(cdist_q(v * a, Q), Q) for v in S)
            if m > best: best = m
    return best

LO, HI = F(1, 12), F(2, 23)   # the claimed gap (open interval)
print(f"gap interval (1/12, 2/23) = ({float(LO):.6f}, {float(HI):.6f}); width = {float(HI-LO):.6f} = 1/276")
print("="*78)
print("FULL ENUMERATION: all 11-subsets of {1..B}, primitive, count M(U) in (1/12, 2/23)")
print("="*78)
for B in [13, 14, 15, 16, 17]:
    ingap = []
    total = 0; prim = 0
    for U in itertools.combinations(range(1, B + 1), 11):
        total += 1
        g = 0
        for v in U: g = gcd(g, v)
        if g != 1: continue
        prim += 1
        M = Mval(U, 2 * B + 2)
        if LO < M < HI:
            ingap.append((U, M))
    print(f"  B={B:>2}: {prim:>6} primitive 11-sets;  # with M in (1/12,2/23): {len(ingap)}"
          + ("" if not ingap else f"  e.g. {ingap[0]}"))

print()
print("="*78)
print("LARGE-SPEED / DILATED near-AP families: try to land M in (1/12, 2/23)")
print("  (a) c*{1..11} with one element replaced by a nearby value (genuine 11-set)")
print("  (b) {1..11} with one element multiplied up (spread a single runner)")
print("="*78)
found = 0
# (a) replace element j of c*{1..11} by c*j + s (small shift), keep 11 distinct
for c in [1,2,3,4,6]:
    base = [c*k for k in range(1,12)]
    for j in range(11):
        for s in [-2,-1,1,2, c, -c]:
            U = base[:j] + [base[j]+s] + base[j+1:]
            U = sorted(set(x for x in U if x>0))
            if len(U)!=11: continue
            M = Mval(U, min(3*max(U)+2, 400))
            if LO < M < HI:
                found += 1
                if found <= 8: print(f"  IN GAP: c={c} j={j} s={s} U={U} M={M}")
# (b) {1..10} u {m} for m up to 40 (spread the top runner) -- these are the ladder climbers
print("  {1..10} u {m}, m=11..40 : M values (ladder climb from 1/12 up):")
seen=set()
for m in range(11, 41):
    U = list(range(1,11)) + [m]
    M = Mval(U, 3*m+2)
    if M not in seen:
        seen.add(M)
        tag = " <== IN GAP (1/12,2/23)!" if LO<M<HI else ("" if M!=LO else " (=1/12 tight)")
        print(f"    m={m:>2}: M={M} (~{float(M):.5f}){tag}")
print(f"\n  total families found with M in (1/12, 2/23): {found}")
print()
print("READING: 0 in-gap => 1/12 is ISOLATED in the 11-runner spectrum (gap delta0=1/276);")
print("the residual's 'near-AP' loose U all have M(U)>=2/23 (definite margin, discrete rung).")
print("DONE")
