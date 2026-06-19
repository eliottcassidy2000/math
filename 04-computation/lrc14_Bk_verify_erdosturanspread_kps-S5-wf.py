#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_erdosturanspread_kps-S5-wf.py   (kind-pasteur-2026-06-18-S5, ADVERSARIAL VERIFIER)

Independent adversarial verification of the "erdos-turan-spread-bound" claimed advance.
We DO NOT import the claimed code; we re-implement mu_theta(E) from scratch with a DIFFERENT
rigorous method, then cross-check every quoted constant and HUNT for counterexamples.

mu_theta(E) = meas{ x in [0,1) : circular maxgap of {frac(e x): e in E} > theta }.
E = sorted nonneg ints, 0 in E, k=|E|, spread = max E.

INDEPENDENT EXACT ENGINE (different from claimed):
  On each order-cell (a,b) between collision breakpoints (e_i-e_j)x in Z, the cyclic order of
  points is fixed and each point p_i(x)=e_i x - floor(e_i*mid) is affine. Each cyclic gap g(x)
  is affine. {maxgap>theta} = UNION over gaps of {g>theta}. We compute that union's measure
  exactly per cell. (Logically identical math, INDEPENDENT implementation, with redundant
  cross-checks: we also verify via a SECOND method -- direct gap-crossing breakpoints -- that
  the per-cell measure agrees.)
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass

# ----------------------------------------------------------------------------
# ENGINE 1: order-cell + union-of-affine-half-intervals (rigorous)
# ----------------------------------------------------------------------------
def collision_breakpoints(E):
    bp = {F(0), F(1)}
    Es = sorted(set(E))
    for i in range(len(Es)):
        for j in range(i+1, len(Es)):
            d = Es[j]-Es[i]
            for m in range(0, d+1):
                bp.add(F(m, d))
    return sorted(b for b in bp if F(0) <= b <= F(1))

def cell_cyclic_gaps(E, a, b):
    """Affine cyclic gaps (alpha,beta) on open cell (a,b): each gap = alpha*x+beta."""
    mid = (a+b)/2
    Es = sorted(set(E))
    pts = []
    for e in Es:
        val = e*mid
        n = val - (val % 1)          # exact floor
        pts.append((val - n, e, n))  # (frac at mid, e, integer level)
    pts.sort(key=lambda t: t[0])
    m = len(pts)
    gaps = []
    for i in range(m):
        (_, e_i, n_i) = pts[i]
        if i < m-1:
            (_, e_j, n_j) = pts[i+1]
            alpha = e_j - e_i
            beta  = n_i - n_j
        else:
            (_, e0, n0) = pts[0]
            alpha = e0 - e_i
            beta  = n_i - n0 + 1
        gaps.append((F(alpha), F(beta)))
    return gaps

def union_gt(gaps, a, b, theta):
    """measure{ x in (a,b): max_g(alpha x+beta) > theta } = union_g{g>theta} cap (a,b)."""
    ivs = []
    for alpha, beta in gaps:
        if alpha == 0:
            if beta > theta: ivs.append((a, b))
        else:
            xstar = (theta - beta)/alpha
            if alpha > 0: lo, hi = max(a, xstar), b
            else:         lo, hi = a, min(b, xstar)
            if lo < hi: ivs.append((lo, hi))
    if not ivs: return F(0)
    ivs.sort()
    tot = F(0); clo, chi = ivs[0]
    for lo, hi in ivs[1:]:
        if lo <= chi: chi = max(chi, hi)
        else: tot += chi-clo; clo, chi = lo, hi
    tot += chi-clo
    return tot

def mu_theta(E, theta):
    E = sorted(set(E))
    if len(E) == 1: return F(1)
    bps = collision_breakpoints(E)
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        tot += union_gt(cell_cyclic_gaps(E, a, b), a, b, theta)
    return tot

# ----------------------------------------------------------------------------
# ENGINE 2 (cross-check): explicit gap=theta breakpoints + midpoint sampling.
# Add, to the collision breakpoints, every x where some affine gap hits theta.
# Then on each refined sub-cell the truth of {maxgap>theta} is constant -> sample midpoint.
# Fully rigorous and independent of ENGINE 1's union logic.
# ----------------------------------------------------------------------------
def mu_theta_v2(E, theta):
    E = sorted(set(E))
    if len(E) == 1: return F(1)
    bps = set(collision_breakpoints(E))
    base = sorted(bps)
    # within each collision-cell, add gap=theta crossings
    extra = set()
    for a, b in zip(base, base[1:]):
        if a == b: continue
        for alpha, beta in cell_cyclic_gaps(E, a, b):
            if alpha != 0:
                xstar = (theta - beta)/alpha
                if a < xstar < b: extra.add(xstar)
    allbp = sorted(bps | extra)
    tot = F(0)
    for a, b in zip(allbp, allbp[1:]):
        if a == b: continue
        mid = (a+b)/2
        pts = sorted((e*mid) % 1 for e in E)
        if len(pts) == 1: mg = F(1)
        else:
            gaps = [pts[t+1]-pts[t] for t in range(len(pts)-1)] + [pts[0]+1-pts[-1]]
            mg = max(gaps)
        if mg > theta: tot += b-a
    return tot

TWO7 = F(2,7); ONE7 = F(1,7)

def gcd1(E): return reduce(gcd, E) == 1

def Fk(k, L):
    s = F(0); j = 1
    while 1 - j*L > 0:
        s += (-1)**(j+1)*comb(k, j)*(1 - j*L)**(k-1); j += 1
    return s

# ============================================================================
print("="*84)
print("CHECK A: two independent engines agree, and match canon consecutive 2/7 values")
print("="*84)
# Canon claims (from prompt): mu_3=1, mu_4=19/21, mu_5=9/14, mu_6=4/7,
#   CORRECTED consecutive mu_7=83/210 (NOT 13/35; 13/35 is perforated minimizer),
#   mu_13(consec 2/7)=829/4620.
canon27 = {3:F(1),4:F(19,21),5:F(9,14),6:F(4,7),7:F(83,210),13:F(829,4620)}
allok = True
for k in range(3, 14):
    E = list(range(k))
    m1 = mu_theta(E, TWO7)
    m2 = mu_theta_v2(E, TWO7)
    agree = (m1 == m2)
    tag = ""
    if k in canon27:
        ok = (m1 == canon27[k]); allok &= ok
        tag = f"  canon={canon27[k]} {'OK' if ok else 'MISMATCH'}"
    print(f"  k={k:2d} consec  mu_2/7={m1}={float(m1):.6f}  engines_agree={agree}{tag}")
    if not agree: allok = False
print(f"  ==> ALL consecutive 2/7 checks (both engines + canon): {'PASS' if allok else 'FAIL'}")
print(f"  NOTE the prompt's CORRECTION: consecutive mu_7 = 83/210 = {float(F(83,210)):.5f}, "
      f"and mu({{0,2,3,4,5,6,8}})={mu_theta([0,2,3,4,5,6,8],TWO7)} (the claimed 13/35 perforated min)")

print()
print("="*84)
print("CHECK B: F(k) iid-ceiling constants quoted in prompt")
print("="*84)
for k,claim in [(4,F(342,343)),(7,F(13443,16807)),(9,F(3279513,5764801)),(13,F(3132376013,13841287201))]:
    got = Fk(k, TWO7)
    print(f"  F({k}) = {got}  quoted={claim}  {'OK' if got==claim else 'MISMATCH'}")

print()
print("="*84)
print("CHECK C: mu_{1/7} = 1 for k<=7 (pigeonhole claim) -- EXHAUSTIVE over all E up to spread")
print("="*84)
# Claim: mu_{1/7}(E)=1 for EVERY E with k<=7. Test exhaustively over a range of spreads.
c_ok = True
for k in range(3, 8):
    worst = F(1); worstE = None
    smax = k + 8  # generous
    for s in range(k-1, smax+1):
        if comb(s-1, max(k-2,0)) > 60000:
            continue
        for interior in itertools.combinations(range(1, s), k-2):
            E = (0,) + interior + (s,)
            m = mu_theta(list(E), ONE7)
            if m < worst: worst = m; worstE = E
    ok = (worst == F(1)); c_ok &= ok
    print(f"  k={k}: min mu_1/7 over spreads<= {smax} = {worst} {'(=1 OK)' if ok else f'<1 at {worstE}'}")
print(f"  ==> mu_1/7==1 for all tested k<=7: {'CONFIRMED' if c_ok else 'REFUTED'}")
