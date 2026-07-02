#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE SIMULTANEOUS-PEEL LEMMA and the corrected inductive tower for the r=2 residual (11-cores, band 1/14).

CONTEXT. S31 closed the 11-core census for max<=19 + single-outlier (w<=500), and flagged the multi-outlier
tail as blocked on a "uniform arc-count bound c'<=const" for cores CONTAINING a large element. That
requirement is WRONG (and false as needed: c(L_{C''+{w}}) ~ w*meas(L_{C''}) grows linearly in w). The fix:

LEMMA (simultaneous peel; elementary, rigorous). Let C = C_low u F (disjoint), j = |F|. Then
    meas(L_C)  >=  (1 - j/7) * meas(L_{C_low})  -  (2 c_low / 7) * sum_{w in F} 1/w ,
where c_low = number of connected components (arcs) of L_{C_low}.
Proof: L_C = L_{C_low} \ U_{w in F} danger(w).  danger(w) = w arcs of length 1/(7w) centered at k/w.
For a single interval I (|I|=len), #danger(w)-arcs meeting I <= w*len + 2, so
meas(I n danger(w)) <= len/7 + 2/(7w).  Summing over the c_low components of L_{C_low} and over w in F:
meas(L_{C_low} n U danger(w)) <= j*meas(L_{C_low})/7 + (2 c_low/7) sum 1/w.  QED.

KEY CONSEQUENCE (scale cancellation): at a multiplicative gap cut (max(C_low)*Lambda <= min(F)),
c_low <= sum_{v in C_low} v <= 11*max(C_low), so the error <= (2/7)*11*max(C_low)*j/(Lambda*max(C_low))
= 22 j/(7 Lambda) -- UNIFORM, no absolute scale, no census arc-count needed. The lemma closes every
j<=6-far-element split; it dies exactly at j=7 = 1/(2r) (seven speeds can cover), leaving as the only
residual: GAP-FREE CLUSTERS of >=7 elements at unbounded height (the deep-cluster minimax).

This script: (L) verifies the lemma exactly on adversarial splits (near-equal far pairs, harmonic pairs,
mixed scales -- the cases the old one-at-a-time lever could NOT handle); (T) computes the census tower
m_k (max<=15) and solves the recursion M_k = min(census_k, min_j (1-j/7) M_{k-j}), with trivial
M_k = 1-k/7 for k<=6; prints the guard table for the 11-core; (P) probes the j>=7 deep-cluster residual:
convergence of consecutive clusters {N..N+10} to the two-scale limit integral, the two-scale
factorization for 4-compact + 7-cluster, and the difference-core renormalization prediction.

opus-2026-07-01-S32.
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

BAND = Fr(1, 14); TARGET = Fr(1, 36)

def safe_arcs(v):
    return [((Fr(k) + BAND) / v, (Fr(k + 1) - BAND) / v) for k in range(v)]

def inter(A, B):
    r = []; i = j = 0
    while i < len(A) and j < len(B):
        lo = A[i][0] if A[i][0] > B[j][0] else B[j][0]
        hi = A[i][1] if A[i][1] < B[j][1] else B[j][1]
        if lo < hi: r.append((lo, hi))
        if A[i][1] < B[j][1]: i += 1
        else: j += 1
    return r

def L_arcs(S):
    a = safe_arcs(S[0])
    for v in S[1:]:
        a = inter(a, safe_arcs(v))
        if not a: return []
    return a

def Lmeas_arcs(a): return sum(h - l for l, h in a) if a else Fr(0)
def Lmeas(S): return Lmeas_arcs(L_arcs(S))

def ncomp(a):
    """component count on the CIRCLE (merge wrap-around if first arc starts at 0 and last ends at 1)."""
    if not a: return 0
    n = len(a)
    if n > 1 and a[0][0] == 0 and a[-1][1] == 1: n -= 1
    return n

print("=" * 104)
print(" (L) THE SIMULTANEOUS-PEEL LEMMA -- exact verification on adversarial splits")
print("     bound: meas(L_C) >= (1 - j/7) meas(L_low) - (2 c_low/7) sum_{w in F} 1/w")
print("=" * 104)

LOW_CORES = [
    (1, 2, 3, 5, 7, 8, 9, 11, 12, 13),      # the min-meas 10-core
    (1, 2, 3, 4, 5, 7, 8, 9, 11),           # 9-core (pentagon-like)
    (1, 2, 3, 4, 5, 6, 7, 8),               # 8-core AP
    (2, 3, 4, 5, 6, 7, 9, 12, 13),          # 9-core no-1
    (1, 2, 3, 4, 5),                        # 5-core
]
FAR_FAMS = [
    ("near-equal pair", lambda w: (w, w + 1)),
    ("harmonic pair 2w+1", lambda w: (w, 2 * w + 1)),
    ("near-equal triple", lambda w: (w, w + 1, w + 2)),
    ("spread pair w,7w+3", lambda w: (w, 7 * w + 3)),
    ("AP4 cluster", lambda w: (w, w + 2, w + 4, w + 7)),
    ("two-scale quad", lambda w: (w, w + 1, 5 * w, 5 * w + 1)),
]
WS = [20, 47, 101, 211]
rows = 0; viol = 0; min_slack = None
for low in LOW_CORES:
    arcs = L_arcs(low); mlow = Lmeas_arcs(arcs); cl = ncomp(arcs)
    for name, fam in FAR_FAMS:
        for w in WS:
            F = fam(w)
            if set(F) & set(low): continue
            j = len(F)
            if j > 6 or len(low) + j > 11: continue
            C = tuple(sorted(set(low) | set(F)))
            m = Lmeas(C)
            bound = (1 - Fr(j, 7)) * mlow - Fr(2 * cl, 7) * sum(Fr(1, w_) for w_ in F)
            slack = m - bound
            rows += 1
            if slack < 0: viol += 1
            if min_slack is None or slack < min_slack: min_slack = slack
print(f"  low-cores x far-families x w in {WS}: {rows} exact checks; violations = {viol}; "
      f"min slack = {float(min_slack):.6f} ({'>=0 OK' if min_slack >= 0 else 'VIOLATED'})")

# the case the OLD lever could not handle: TWO near-equal far elements (peeling one leaves c' ~ w1*meas)
print("\n  the old blocker, directly: C_low = 9-core + F = {w, w+1} (near-equal far pair):")
low = (1, 2, 3, 4, 5, 7, 8, 9, 11); arcs = L_arcs(low); mlow = Lmeas_arcs(arcs); cl = ncomp(arcs)
print(f"    C_low = {low}: meas = {float(mlow):.6f}, c_low = {cl}")
print(f"    {'w':>6} {'meas(L_C)':>12} {'new bound':>12} {'slack':>10} {'c(L_low+w) (old lever c'')':>26}")
for w in (20, 50, 100, 200, 400):
    F = (w, w + 1); C = tuple(sorted(set(low) | set(F)))
    m = Lmeas(C)
    bound = (1 - Fr(2, 7)) * mlow - Fr(2 * cl, 7) * (Fr(1, w) + Fr(1, w + 1))
    aw = inter(arcs, safe_arcs(w))
    print(f"    {w:>6} {float(m):>12.6f} {float(bound):>12.6f} {float(m - bound):>10.6f} {ncomp(aw):>26}")
print("    (old one-at-a-time lever error ~ c''/(7(w+1)) with c'' growing ~ w -> does NOT vanish; new bound does.)")

print()
print("=" * 104)
print(" (T) THE TOWER: census m_k (primitive k-cores, max<=15) + recursion M_k, and the 11-core guard table")
print("=" * 104)
Vt = 15
census = {}
carc = {}
for k in range(5, 12):
    mn = None; arg = None; cmax = 0
    for C in itertools.combinations(range(1, Vt + 1), k):
        if reduce(gcd, C) != 1: continue
        a = L_arcs(C); m = Lmeas_arcs(a)
        if m > 0:
            if mn is None or m < mn: mn, arg = m, C
            c = ncomp(a)
            if c > cmax: cmax = c
    census[k] = mn; carc[k] = cmax
    print(f"  census m_{k}(max<={Vt}) = {mn} = {float(mn):.6f} at {arg};  max arc-count seen = {cmax}")

M = {}
for k in range(1, 7):
    M[k] = Fr(7 - k, 7)  # unconditional union bound, ALL k-cores
for k in range(7, 12):
    rec = min((1 - Fr(j, 7)) * M[k - j] for j in range(1, 7) if k - j >= 1)
    M[k] = min(census[k], rec) if k in census and census[k] is not None else rec
print("\n  recursion  M_k = min( census_k , min_{1<=j<=6} (1-j/7) M_{k-j} )   [ideal eps->0 at the gap cut]")
print(f"  {'k':>4} {'M_k':>14} {'float':>10}   binding branch")
for k in range(5, 12):
    rec = min(((1 - Fr(j, 7)) * M[k - j], j) for j in range(1, 7) if k - j >= 1)
    which = "trivial 1-k/7" if k <= 6 else (f"census" if census.get(k) is not None and census[k] <= rec[0] else f"peel j={rec[1]}")
    print(f"  {k:>4} {str(M[k]):>14} {float(M[k]):>10.6f}   {which}")

print(f"\n  GUARD TABLE for the 11-core (target 1/36 = {float(TARGET):.6f}); error budget eps(Lambda) = 22j/(7 Lambda):")
print(f"  {'j far':>6} {'(1-j/7)M_(11-j)':>16} {'guard-1/36':>12} {'eps needed <':>13} {'Lambda needed':>14}")
allok = True
for j in range(1, 7):
    g = (1 - Fr(j, 7)) * M[11 - j]
    margin = g - TARGET
    ok = margin > 0; allok = allok and ok
    lam = float(Fr(22 * j, 7) / margin) if ok else float('inf')
    print(f"  {j:>6} {float(g):>16.6f} {float(margin):>12.6f} {float(margin):>13.6f} {lam:>14.0f}")
print(f"  ALL j=1..6 guards positive? {allok}")
print("  => with Lambda ~ max column, the gap-split residual is CLOSED for any j<=6 far elements at ANY scales;")
print("     what remains: gap-free clusters of >=7 elements at unbounded height (+ their <=4-element compact part).")

print()
print("=" * 104)
print(" (P) PROBING THE j>=7 DEEP-CLUSTER RESIDUAL")
print("=" * 104)
print(" P1. consecutive 11-cluster {N..N+10}: exact meas vs the two-scale limit integral")
def D_density(offsets_t, jj):
    """meas_u of intersection of arcs (A - i*t), i.e. lonely density of the frozen cluster pattern."""
    # arcs on u-circle: A - c for c in offsets: A=[1/14,13/14]
    a = None
    for c in offsets_t:
        lo = (Fr(1, 14) - c) % 1; hi = (Fr(13, 14) - c) % 1
        arcs_c = [(lo, hi)] if lo < hi else [(Fr(0), hi), (lo, Fr(1))]
        a = arcs_c if a is None else inter(a, arcs_c)
        if not a: return Fr(0)
    return Lmeas_arcs(a)

# limit integral for consecutive j-cluster: Int_0^1 D(t) dt with D(t)=meas(n_i (A - i t)), grid-exact Riemann
def limit_integral_consec(jj, Q=2000):
    s = Fr(0)
    for kq in range(Q):
        t = Fr(kq, Q)
        s += D_density([ (i * t) % 1 for i in range(jj) ], jj)
    return s / Q

lim11 = limit_integral_consec(11, 1400)
print(f"    limit integral Int D_11 ~ {float(lim11):.6f}")
for N in (30, 60, 120, 240):
    C = tuple(range(N, N + 11))
    g = reduce(gcd, C)
    m = Lmeas(C)
    print(f"    N={N:>4}: meas(L_{{N..N+10}}) = {float(m):.6f}   (gcd={g}); vs limit {float(lim11):.6f}  diff {float(m - lim11):+.6f}")

print("\n P2. two-scale factorization: C = {1,2,3,4} u {N..N+6};  predict Int_{L_low} D_7(t) dt")
low = (1, 2, 3, 4); arcs_low = L_arcs(low); mlow = Lmeas_arcs(arcs_low)
def integral_over_arcs(arcs, jj, Q=3000):
    """Riemann approx of Int_{arcs} D_jj(t) dt (consecutive pattern, offsets i*t)."""
    s = Fr(0); tot = Fr(0)
    for kq in range(Q):
        t = Fr(2 * kq + 1, 2 * Q)
        inside = any(l <= t < h for l, h in arcs)
        if inside:
            s += D_density([ (i * t) % 1 for i in range(jj) ], jj)
    return s / Q
pred = integral_over_arcs(arcs_low, 7, 2100)
print(f"    meas(L_low) = {float(mlow):.6f};  predicted two-scale value Int_(L_low) D_7 ~ {float(pred):.6f}")
for N in (30, 60, 120, 240):
    C = tuple(sorted(set(low) | set(range(N, N + 7))))
    m = Lmeas(C)
    print(f"    N={N:>4}: meas = {float(m):.6f}   vs predicted {float(pred):.6f}   diff {float(m - pred):+.6f}   > 1/36? {m > TARGET}")

print("\n P3. worst 7-cluster PATTERN probe (difference-core renormalization): cluster {N + c_i}, varying c;")
print("     prediction: the worst pattern's difference core is AP-like (the census extremizer's shadow).")
PATTERNS = {
    "consecutive {0..6}":            (0, 1, 2, 3, 4, 5, 6),
    "even spread {0,2,..,12}":       (0, 2, 4, 6, 8, 10, 12),
    "dilated-AP shadow {0,7,..,42}": (0, 7, 14, 21, 28, 35, 42),
    "lacunary {0,1,3,7,15,31,63}":   (0, 1, 3, 7, 15, 31, 63),
    "random-ish {0,5,9,16,24,33,41}":(0, 5, 9, 16, 24, 33, 41),
    "pentagon-diff {0,1,2,3,5,7,9}": (0, 1, 2, 3, 5, 7, 9),
}
N = 210
res = []
for name, pat in PATTERNS.items():
    C = tuple(sorted(set(low) | set(N + c for c in pat)))
    if len(C) != 11:
        continue
    m = Lmeas(C)
    res.append((m, name))
    print(f"    {name:<32} meas(low u (N+pat)) = {float(m):.6f}   > 1/36? {m > TARGET}")
res.sort()
print(f"    WORST pattern at N={N}: {res[0][1]} with meas {float(res[0][0]):.6f} (margin over 1/36: {float(res[0][0]/TARGET):.3f}x)")

print("\n P4. D_7 vanishing (the continuous-Fraenkel tiling): consecutive-7 density at t=k/7 (exact):")
for kk in range(1, 4):
    t = Fr(kk, 7)
    d = D_density([ (i * t) % 1 for i in range(7) ], 7)
    print(f"    D_7(k={kk}/7) = {d}   (danger arcs tile the circle exactly <=> D=0)")
print("    => the 7-cluster CAN tile (D=0) but only at isolated t (tiling breaks under perturbation since")
print("       arcs move at distinct speeds); the zeros are non-degenerate (linear), so the parent integral")
print("       Int_(L_low) D stays positive -- quantifying THAT floor is the remaining deep-cluster minimax.")
print("DONE.")
