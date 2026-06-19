#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_unionfloor_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL)

Verify the k>=8 union-bound floor INDEPENDENTLY:
  rho*_{1/7}(P,E) := meas{ x in G_P : maxgap{frac(e x)} > 1/7 } >= meas(G_P) + mu_{1/7}(E) - 1.
This is just inclusion-exclusion (both subsets of [0,1)): meas(A cap B) >= meas A + meas B - 1.
=> if for EVERY admissible (P,E) we have meas(G_P) + mu_{1/7}(E) > 1, then rho*_{1/7} > 0.

We compute, k=8..13:
  min over admissible P (|P|=13-k) of meas(G_P), call it gp_min(k);
  GRANT the 1/7 spread bound mu_{1/7}(E) >= mu_{1/7}(consec_k) [verified separately];
  the per-k union floor = gp_min(k) + consec_mu17(k) - 1.
We check this is > 0 for all k>=8 and identify the binding k, and verify the quoted
union floor 1891/5880 and threshold thr_8=3637/5880.

CAUTION on what union bound actually delivers: it requires meas(G_P) + mu_{1/7}(E) - 1 > 0
for the SPECIFIC E that arises, with the SPECIFIC P. We test min_P over ALL admissible P, and
GRANT mu_{1/7}(E)>=consec. The binding case is the SMALLEST meas(G_P) paired with consec_mu.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass

H = F(1,14); ONE7 = F(1,7)

def merge(iv):
    iv = sorted(iv); out = []
    for a,b in iv:
        if out and a <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], b))
        else: out.append((a,b))
    return out
def meas(arcs): return sum((b-a for a,b in arcs), F(0))
def complement(arcs):
    arcs = merge(arcs); out = []; prev = F(0)
    for a,b in arcs:
        if a > prev: out.append((prev,a))
        prev = max(prev,b)
    if prev < 1: out.append((prev, F(1)))
    return out
def danger_arcs(u, h=H):
    """{x: ||u x|| < h} = union over residues j of (j/u +- h/u)."""
    iv = []
    for j in range(u):
        c = F(j,u); a = (c - h/u) % 1; b = (c + h/u) % 1
        if a < b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def safe_set(P, h=H):
    """G_P = {x: ||p x|| >= 1/14 for all p in P}."""
    if not P: return [(F(0),F(1))]
    return complement(merge([iv for u in P for iv in danger_arcs(u,h)]))

# independent mu_1/7 engine (reuse the verified one)
def collision_breakpoints(E):
    bp = {F(0), F(1)}; Es = sorted(set(E))
    for i in range(len(Es)):
        for j in range(i+1, len(Es)):
            d = Es[j]-Es[i]
            for m in range(0, d+1): bp.add(F(m, d))
    return sorted(b for b in bp if F(0) <= b <= F(1))
def cell_cyclic_gaps(E, a, b):
    mid = (a+b)/2; Es = sorted(set(E)); pts = []
    for e in Es:
        val = e*mid; n = val - (val % 1); pts.append((val - n, e, n))
    pts.sort(key=lambda t: t[0]); m = len(pts); gaps = []
    for i in range(m):
        (_, e_i, n_i) = pts[i]
        if i < m-1: (_, e_j, n_j) = pts[i+1]; alpha = e_j-e_i; beta = n_i-n_j
        else: (_, e0, n0) = pts[0]; alpha = e0-e_i; beta = n_i-n0+1
        gaps.append((F(alpha), F(beta)))
    return gaps
def union_gt(gaps, a, b, theta):
    ivs = []
    for alpha, beta in gaps:
        if alpha == 0:
            if beta > theta: ivs.append((a,b))
        else:
            xstar = (theta-beta)/alpha
            if alpha > 0: lo, hi = max(a, xstar), b
            else: lo, hi = a, min(b, xstar)
            if lo < hi: ivs.append((lo,hi))
    if not ivs: return F(0)
    ivs.sort(); tot = F(0); clo, chi = ivs[0]
    for lo, hi in ivs[1:]:
        if lo <= chi: chi = max(chi, hi)
        else: tot += chi-clo; clo, chi = lo, hi
    tot += chi-clo; return tot
def mu17(E):
    E = sorted(set(E))
    if len(E) == 1: return F(1)
    bps = collision_breakpoints(E); tot = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        tot += union_gt(cell_cyclic_gaps(E, a, b), a, b, ONE7)
    return tot

print("="*88)
print("CHECK: gp_min(k) = min over admissible P (|P|=13-k, P subset {1..13}) of meas(G_P)")
print("="*88)
gp_min = {}
for k in range(8, 14):
    psz = 13 - k
    if psz == 0:
        gp_min[k] = F(1); print(f"  k={k:2d}: |P|=0 -> meas(G_P)=1"); continue
    best = None; bestP = None
    for P in itertools.combinations(range(1, 14), psz):
        m = meas(safe_set(list(P)))
        if best is None or m < best: best = m; bestP = P
    gp_min[k] = best
    print(f"  k={k:2d}: |P|={psz}  gp_min = {best} = {float(best):.6f}  at P={bestP}")

print()
print("="*88)
print("CHECK: per-k union floor = gp_min(k) + mu_1/7(consec_k) - 1  (GRANTING the spread bound)")
print("="*88)
consec = {k: mu17(list(range(k))) for k in range(8,14)}
floors = {}
for k in range(8, 14):
    fl = gp_min[k] + consec[k] - 1
    floors[k] = fl
    print(f"  k={k:2d}: gp_min={float(gp_min[k]):.5f} + consec_mu17={float(consec[k]):.5f} - 1 = {fl} = {float(fl):.6f}  {'>0 OK' if fl>0 else '<=0 FAIL'}")
gmin = min(floors.values())
binding = [k for k in floors if floors[k]==gmin]
print(f"\n  GLOBAL k>=8 union floor = {gmin} = {float(gmin):.6f}  (binding k={binding})")
print(f"  quoted floor 1891/5880 = {F(1891,5880)} = {float(F(1891,5880)):.6f}  "
      f"{'MATCH' if gmin==F(1891,5880) else 'DIFFERS'}")

print()
print("="*88)
print("CHECK: quoted thr_8 = 3637/5880 = 1 - gp_min(8)?  (thr_k := 1 - min meas(G_P))")
print("="*88)
thr8 = 1 - gp_min[8]
print(f"  1 - gp_min(8) = {thr8} = {float(thr8):.6f}   quoted thr_8=3637/5880={float(F(3637,5880)):.6f}  "
      f"{'MATCH' if thr8==F(3637,5880) else 'DIFFERS'}")
print(f"  consec_mu17(8) = {consec[8]} = {float(consec[8]):.5f} ; margin consec - thr_8 = {consec[8]-thr8} = {float(consec[8]-thr8):.5f}")
print(f"  (quoted margin 1891/5880 = {float(F(1891,5880)):.5f})")

print()
print("="*88)
print("ADVERSARIAL: is the union bound's hypothesis met by the WORST admissible (P,E) pair?")
print("The union bound needs meas(G_P)+mu_1/7(E)-1>0 for the actual E (any E, not just consec).")
print("GRANTING mu_1/7(E)>=consec, worst pair = (P*=argmin G_P, E=consec). Already covered above.")
print("But cross-check: directly compute rho*_1/7(P*,consec_E) and confirm >= union floor.")
print("="*88)
def good17_arcs(E):
    E = sorted(set(E))
    bps = collision_breakpoints(E); good = []
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        m = union_gt(cell_cyclic_gaps(E, a, b), a, b, ONE7)
        # need the actual arcs, not just measure -> recompute as intervals
        gaps = cell_cyclic_gaps(E, a, b); ivs = []
        for alpha, beta in gaps:
            if alpha == 0:
                if beta > ONE7: ivs.append((a,b))
            else:
                xstar = (ONE7-beta)/alpha
                if alpha > 0: lo, hi = max(a,xstar), b
                else: lo, hi = a, min(b,xstar)
                if lo < hi: ivs.append((lo,hi))
        good.extend(ivs)
    return merge(good)
def intersect(A, B):
    A = merge(A); B = merge(B); out = []; i=j=0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi: out.append((lo,hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return out
for k in range(8, 14):
    psz = 13-k
    # find argmin P
    bestP = None; best = None
    for P in itertools.combinations(range(1,14), psz):
        m = meas(safe_set(list(P)))
        if best is None or m < best: best = m; bestP = P
    GP = safe_set(list(bestP)) if bestP else [(F(0),F(1))]
    GE = good17_arcs(list(range(k)))
    rho = meas(intersect(GP, GE))
    ub = floors[k]
    print(f"  k={k:2d}: actual rho*_1/7(P*,consec) = {rho} = {float(rho):.6f}  ;  union LB = {float(ub):.6f}  "
          f"{'(rho>=LB OK)' if rho>=ub else '(rho<LB!!)'}")
print("\nDONE.")
