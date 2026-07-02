#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3950: the MULTI-OUTLIER residual of the r=2 moment relaxation (11-core floor inf meas >= 1/36).

kind-pasteur-2026-07-01-S28. State: kps-S27 census closed the BOUNDED 11-core case (min = pentagon
313/9702 = 0.032261, all >= 1/36, margin 1.161x). Single-outlier (10 bounded + 1 far) min = 0.03231
at {1,2,3,5,7,8,9,11,12,13}+50 = the (6/7) x (10-core min) equidistribution limit. THE RESIDUAL is
the multi-outlier case (>= 2 far speeds). This script:

 (1) BASE LEDGER: min over bounded m-cores of meas(L_B), m = 7..11, with the arc count N_B of each
     minimizer. Tests the FLAT-RATIO LAW min_m ~ (7/6)^(11-m) * 0.0323 (per-slot cost = one
     equidistribution factor) and the GREEDY-DROP NESTING of the argmins.
 (2) THE UNIFORM ARC-COUNT UNION FLOOR (the elementary theorem this session contributes):
         meas(L_B ∩ ∩_j safe(w_j)) >= meas(L_B)*(1 - r/7) - (N_B/7) * Σ_j 1/w_j
     valid for ALL integer w_j (periodicity repairs the union bound: a 1/w-periodic comb cannot
     concentrate in L_B beyond meas/7 + N_B/(7w)). Explicit finite cutoffs W*(r) follow.
 (3) ONE-STEP PEEL + ARC GROWTH verification:
         meas(L ∩ safe_w) >= (6/7) meas(L) - N(L)/(7w),   N(L∩safe_w) <= 2N(L) + w*meas(L)
 (4) SINGLE-OUTLIER scan (verify 0.03231 @ {1..13}\{4,6,10}+50), 2-OUTLIER and 3-OUTLIER scans
     including resonant structures (w2=2w1, w2=w1+d, arithmetic triples), vs the flat prediction.
 (5) PAIR OVERLAP: integer combs cannot be disjoint; meas(D_w ∩ D_w') vs 1/49 (independent) and the
     resonant excess (e.g. (w,2w): overlap = 1/14).
Verdict: where union floor + ledger closes (r<=?), where the peel/spread machinery is needed, and
whether any multi-outlier configuration dips below the single-outlier min or 1/36.
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import gcd, floor
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

BAND = 1.0/14.0
T36  = 1.0/36.0

# ---------------------------------------------------------------- float engine (per-arc comb clip)
def clip(arcs, v, band=BAND):
    """intersect arc list with safe(v) = {t : ||vt|| > band} (teeth [ (k+band)/v, (k+1-band)/v ])."""
    out = []
    for (lo, hi) in arcs:
        k0 = int(floor(lo*v)); k1 = int(floor(hi*v))
        for k in range(k0, k1+1):
            a = (k + band)/v; b = (k + 1 - band)/v
            l = lo if lo > a else a
            h = hi if hi < b else b
            if l < h - 1e-15: out.append((l, h))
    return out

def lonely_arcs(S, band=BAND):
    arcs = [(0.0, 1.0)]
    for v in sorted(S):
        arcs = clip(arcs, v, band)
        if not arcs: return arcs
    return arcs

def measN(S):
    a = lonely_arcs(S)
    return sum(h-l for l, h in a), len(a)

# ---------------------------------------------------------------- exact engine (verification)
def measL_exact(S):
    BQ = Fr(1, 14)
    S = sorted(set(S))
    def safe(v): return [((Fr(k)+BQ)/v, (Fr(k+1)-BQ)/v) for k in range(v)]
    a = safe(S[0])
    for v in S[1:]:
        r = []; i = j = 0; b = safe(v)
        while i < len(a) and j < len(b):
            lo = a[i][0] if a[i][0] > b[j][0] else b[j][0]
            hi = a[i][1] if a[i][1] < b[j][1] else b[j][1]
            if lo < hi: r.append((lo, hi))
            if a[i][1] < b[j][1]: i += 1
            else: j += 1
        a = r
        if not a: return Fr(0)
    return sum(h-l for l, h in a)

def primitive(C):
    g = 0
    for c in C: g = gcd(g, c)
    return tuple(sorted(c//g for c in C))

# =================================================================================================
print("="*104)
print(" (1) BASE LEDGER: min over BOUNDED m-cores, m = 7..11  (drops of {1..V} + random primitive)")
print("="*104)
rng = random.Random(28)
ledger = {}
for m in range(7, 12):
    pool = set()
    for V in range(m, min(m+5, 18)):
        d = V - m
        if d == 0: pool.add(tuple(range(1, V+1))); continue
        from math import comb
        if comb(V, d) <= 60000:
            for drop in itertools.combinations(range(1, V+1), d):
                pool.add(primitive(tuple(sorted(set(range(1, V+1)) - set(drop)))))
    for _ in range(4000):
        C = tuple(sorted(rng.sample(range(1, 31), m)))
        pool.add(primitive(C))
    best = []
    for C in pool:
        mu, N = measN(C)
        best.append((mu, N, C))
    best.sort()
    ledger[m] = best
    mu, N, C = best[0]
    r = 11 - m
    flat  = mu * (6.0/7.0)**r
    union = mu * (1.0 - r/7.0)
    print(f"  m={m:2d}  #cores={len(pool):6d}  min meas={mu:.6f}  N_arcs={N:3d}  argmin={C}")
    print(f"         flat (6/7)^{r} * min = {flat:.6f}   union (1-{r}/7)*min = {union:.6f}   vs 1/36={T36:.6f}"
          f"   flat_clears={flat>=T36}  union_clears={union>=T36}")
print("\n  FLAT-RATIO TEST: min_m * (7/6)^(11-m) (should be ~constant if per-slot cost = 6/7):")
for m in range(7, 12):
    mu = ledger[m][0][0]
    print(f"    m={m:2d}: min={mu:.6f}   normalized min*(6/7)^(11-m) = {mu*(6.0/7.0)**(11-m):.6f}")
print("\n  GREEDY-DROP NESTING of argmins (is argmin_{m} = argmin_{m+1} minus one element?):")
chain = [set(ledger[m][0][2]) for m in range(11, 6, -1)]
for i in range(len(chain)-1):
    a, b = chain[i], chain[i+1]
    print(f"    argmin_{11-i} ⊃ argmin_{10-i}? {b.issubset(a)}   dropped: {sorted(a-b) if b.issubset(a) else '(not nested)'}")
print("\n  top-5 per m (near-minimal locus):")
for m in range(7, 12):
    for mu, N, C in ledger[m][:5]:
        print(f"    m={m:2d}  {mu:.6f}  N={N:3d}  {C}")

# =================================================================================================
print("\n" + "="*104)
print(" (2) UNIFORM ARC-COUNT UNION FLOOR + cutoffs:  meas >= meas(L_B)(1-r/7) - (N_B/7) Σ 1/w_j")
print("="*104)
print("  Per-speed lemma (periodicity repairs the union bound):  meas(L_B ∩ D_w) <= meas(L_B)/7 + N_B/(7w).")
print("  PROOF: D_w = w teeth of width 1/(7w), 1/w-periodic; an arc I meets <= |I|w+1 teeth, each")
print("  contributing <= 1/(7w) => danger in I <= |I|/7 + 1/(7w); sum over the N_B arcs. QED (elementary).")
print("  => cutoff W*(r,m): smallest W with min_m*(1-r/7) - (N/7)(r/W) >= 1/36 for ALL outliers >= W:")
for r in range(1, 5):
    m = 11 - r
    mu, N, C = ledger[m][0]
    margin = mu*(1.0 - r/7.0) - T36
    # UNIFORM cutoff: sup over ALL bounded bases B of the per-base cutoff W*_B (bases whose union
    # margin <= 0 cannot be closed by the floor at any W — count them; they need the peel instead)
    Wmax = 0.0; nfail = 0; failex = None
    for muB, NB, CB in ledger[m]:
        mB = muB*(1.0 - r/7.0) - T36
        if mB <= 0:
            nfail += 1
            if failex is None: failex = (muB, NB, CB)
            continue
        WB = (NB/7.0)*r/mB
        if WB > Wmax: Wmax = WB
    if margin > 0:
        Wstar = (N/7.0)*r/margin
        print(f"    r={r}: base min_m={mu:.6f} (m={m}, N={N})  union margin={margin:+.6f}  W*(argmin)={Wstar:8.1f}"
              f"   UNIFORM W*max={Wmax:9.1f}   #bases with margin<=0: {nfail}")
        if failex: print(f"         first failing base: meas={failex[0]:.6f} N={failex[1]} {failex[2]}")
    else:
        print(f"    r={r}: base min_m={mu:.6f} (m={m}, N={N})  union margin={margin:+.6f}  =>  union floor FAILS at the"
              f" minimizer; #bases failing: {nfail}  -> those need the multiplicative peel (6/7)^r / pair credit")

# verify the per-speed lemma numerically on ledger minimizers
print("\n  numeric check of  meas(L_B ∩ D_w) <= meas(L_B)/7 + N_B/(7w)  (worst ratio should be <= 1):")
worst = 0.0; worst_case = None
for m in (9, 10, 11):
    mu, N, C = ledger[m][0]
    LB = lonely_arcs(C)
    for w in list(range(15, 121, 7)) + [50, 137, 250, 499, 1000]:
        Ls = clip(LB, w)
        danger = mu - sum(h-l for l, h in Ls)
        bound  = mu/7.0 + N/(7.0*w)
        ratio = danger/bound if bound > 0 else 0.0
        if ratio > worst: worst, worst_case = ratio, (C, w, danger, bound)
print(f"    worst danger/bound = {worst:.4f}  at core={worst_case[0]} w={worst_case[1]}"
      f"  (danger={worst_case[2]:.6f} bound={worst_case[3]:.6f})  LEMMA {'HOLDS' if worst<=1.0 else 'VIOLATED'}")

# =================================================================================================
print("\n" + "="*104)
print(" (3) ONE-STEP PEEL + ARC GROWTH:  meas(L∩safe_w) >= (6/7)meas(L) - N/(7w);  N' <= 2N + w*meas")
print("="*104)
worstP = -1e9; worstG = 0.0; caseP = caseG = None
for m in (9, 10, 11):
    mu, N, C = ledger[m][0]
    LB = lonely_arcs(C); NB = len(LB)
    for w in list(range(15, 200, 11)) + [50, 100, 313, 997]:
        Ls = clip(LB, w); mu2 = sum(h-l for l, h in Ls); N2 = len(Ls)
        deficit = (6.0/7.0)*mu - NB/(7.0*w) - mu2      # should be <= 0 (bound below actual)
        growth  = N2 / (2.0*NB + w*mu)                 # should be <= 1
        if deficit > worstP: worstP, caseP = deficit, (C, w)
        if growth  > worstG: worstG, caseG = growth, (C, w)
print(f"    peel lower bound: worst (bound - actual) = {worstP:+.2e} at {caseP}  {'HOLDS' if worstP<=1e-12 else 'VIOLATED'}")
print(f"    arc growth:       worst N'/(2N + w*meas) = {worstG:.4f}  at {caseG}  {'HOLDS' if worstG<=1.0 else 'VIOLATED'}")

# =================================================================================================
print("\n" + "="*104)
print(" (4a) SINGLE-OUTLIER scan: B = top minimal 10-cores, w = 15..800  (verify 0.03231 @ +50)")
print("="*104)
best1 = []
for mu, N, B in ledger[10][:5]:
    LB = lonely_arcs(B)
    lim = (6.0/7.0)*mu
    lo = None
    for w in range(max(B)+1, 801):
        if w in B: continue
        Ls = clip(LB, w); m2 = sum(h-l for l, h in Ls)
        best1.append((m2, B, w))
        if lo is None or m2 < lo[0]: lo = (m2, w)
    print(f"  B={B} meas_B={mu:.6f} N_B={N}  limit (6/7)meas={lim:.6f}  min over w<=800: {lo[0]:.6f} at w={lo[1]}")
best1.sort()
print("\n  global single-outlier top-5:")
for m2, B, w in best1[:5]:
    print(f"    meas={m2:.6f}  B={B} + w={w}")
m2, B, w = best1[0]
ex = measL_exact(tuple(B)+(w,))
print(f"  EXACT check of global min: meas({B}+{w}) = {ex} = {float(ex):.6f}   >= 1/36? {ex >= Fr(1,36)}"
      f"   margin = {float(ex)/T36:.4f}x")

# =================================================================================================
print("\n" + "="*104)
print(" (4b) TWO-OUTLIER scan: B = top minimal 9-cores; resonant + grid + spread + random pairs")
print("="*104)
best2 = []
bases9 = ledger[9][:6]
for mu, N, B in bases9:
    LB = lonely_arcs(B)
    flat2 = mu*(6.0/7.0)**2
    lo = None
    def try_pair(w1, w2):
        global lo
        if w1 == w2 or w1 in B or w2 in B: return
        Ls = clip(clip(LB, w1), w2); m2 = sum(h-l for l, h in Ls)
        best2.append((m2, B, w1, w2))
        if lo is None or m2 < lo[0]: lo = (m2, w1, w2)
    lo = None
    # resonant lines
    for w1 in range(max(B)+1, 301):
        try_pair(w1, 2*w1); try_pair(w1, 3*w1)
        for d in range(1, 15): try_pair(w1, w1+d)
    # grid
    for w1 in range(max(B)+1, 81):
        for w2 in range(w1+1, min(241, 4*w1)):
            try_pair(w1, w2)
    # spread
    for w1 in range(max(B)+1, 101, 5):
        for K in (10, 50):
            for j in range(7): try_pair(w1, K*w1+j)
    # random
    for _ in range(3000):
        try_pair(rng.randint(14, 999), rng.randint(14, 999))
    print(f"  B={B} meas_B={mu:.6f}  flat (6/7)^2*meas={flat2:.6f}  min found: {lo[0]:.6f} at (w1,w2)=({lo[1]},{lo[2]})")
best2.sort()
print("\n  global two-outlier top-8:")
for m2, B, w1, w2 in best2[:8]:
    print(f"    meas={m2:.6f}  B={B} + ({w1},{w2})   [w2/w1={w2/w1:.3f}, w2-w1={w2-w1}]")
m2, B, w1, w2 = best2[0]
ex = measL_exact(tuple(B)+(w1, w2))
print(f"  EXACT check of global min: {ex} = {float(ex):.6f}  >= 1/36? {ex >= Fr(1,36)}   margin = {float(ex)/T36:.4f}x")
below = [b for b in best2 if b[0] < T36]
print(f"  #pairs below 1/36: {len(below)}  (out of {len(best2)})")

# =================================================================================================
print("\n" + "="*104)
print(" (4c) THREE-OUTLIER probe: B = top minimal 8-cores; structured + random triples")
print("="*104)
best3 = []
for mu, N, B in ledger[8][:3]:
    LB = lonely_arcs(B)
    flat3 = mu*(6.0/7.0)**3
    lo = None
    def try_tri(ws):
        global lo
        ws = tuple(sorted(set(ws)))
        if len(ws) != 3 or any(w in B for w in ws): return
        Ls = LB
        for w in ws: Ls = clip(Ls, w)
        m3 = sum(h-l for l, h in Ls)
        best3.append((m3, B, ws))
        if lo is None or m3 < lo[0]: lo = (m3, ws)
    lo = None
    for w in range(max(B)+1, 201):
        try_tri((w, 2*w, 3*w)); try_tri((w, 2*w, 4*w)); try_tri((w, w+1, w+2))
        try_tri((w, w+1, 2*w)); try_tri((w, w+13, w+26)); try_tri((w, 2*w+1, 4*w+3))
    for w in range(max(B)+1, 101, 3):
        try_tri((w, 10*w, 100*w)); try_tri((w, 10*w+1, 100*w+7))
    for _ in range(6000):
        try_tri((rng.randint(14, 600), rng.randint(14, 600), rng.randint(14, 600)))
    print(f"  B={B} meas_B={mu:.6f}  flat (6/7)^3*meas={flat3:.6f}  min found: {lo[0]:.6f} at {lo[1]}")
best3.sort()
print("\n  global three-outlier top-6:")
for m3, B, ws in best3[:6]:
    print(f"    meas={m3:.6f}  B={B} + {ws}")
m3, B, ws = best3[0]
ex = measL_exact(tuple(B)+ws)
print(f"  EXACT check of global min: {ex} = {float(ex):.6f}  >= 1/36? {ex >= Fr(1,36)}   margin = {float(ex)/T36:.4f}x")
below3 = [b for b in best3 if b[0] < T36]
print(f"  #triples below 1/36: {len(below3)}  (out of {len(best3)})")

# =================================================================================================
print("\n" + "="*104)
print(" (5) PAIR OVERLAP: integer combs cannot be disjoint. meas(D_w1 ∩ D_w2) vs independent 1/49")
print("="*104)
full = [(0.0, 1.0)]
rows = []
for (w1, w2) in [(20,21),(20,40),(20,60),(20,27),(50,51),(50,100),(50,150),(50,63),(97,194),(97,98),(35,36),(35,70),(101,103),(60,84)]:
    S1 = clip(full, w1); m1 = sum(h-l for l, h in S1)
    S12 = clip(S1, w2);  m12 = sum(h-l for l, h in S12)
    # meas(D1 ∩ D2) = 1 - m(S1) - m(S2) + m(S1∩S2) = 1 - 6/7 - 6/7 + m12
    ov = 1.0 - 2.0*(6.0/7.0) + m12
    rows.append((w1, w2, ov))
    tag = "RESONANT" if abs(ov - 1.0/14.0) < 1e-9 else ("~indep" if abs(ov - 1.0/49.0) < 0.004 else "")
    print(f"    (w1,w2)=({w1:3d},{w2:3d})  overlap={ov:.6f}   [1/49={1/49:.6f}, 1/14={1/14:.6f}]  {tag}")
minov = min(r[2] for r in rows)
print(f"    min overlap over sample = {minov:.6f}  (REAL speeds could achieve 0; integer combs stay near >= ~1/49*(1-eps))")

# =================================================================================================
print("\n" + "="*104)
print(" (6) VERDICT / ASSEMBLED SKELETON")
print("="*104)
gmin = min(best1[0][0], best2[0][0], best3[0][0])
print(f"  observed minima: bounded 0.032261 | 1-outlier {best1[0][0]:.6f} | 2-outlier {best2[0][0]:.6f} | 3-outlier {best3[0][0]:.6f}")
print(f"  ALL >= 1/36 = {T36:.6f}?  {gmin >= T36}")
print("""  SKELETON: inf over 11-cores = min over r=#outliers of [bounded (11-r)-core ledger] x [outlier cost]:
   - r=0: S27 census (min = pentagon 313/9702).
   - each r: outliers > W*(r): UNION FLOOR (uniform arc-count, sec 2) — PROVED elementary;
             outliers <= W*(r): FINITE box — scanned here (sec 4), exhaustible in principle;
             thin-margin r (union floor fails): multiplicative PEEL (sec 3) with spread/cluster split,
             clustered pairs gain the overlap credit (sec 5: integer combs overlap >= ~1/49).
  The residual analytic lemma = the peel error control for CLUSTERED outlier packs (the same
  Dedekind/pair-discrepancy object as THM-563/564) — everything else is finite.""")
print("DONE.")
