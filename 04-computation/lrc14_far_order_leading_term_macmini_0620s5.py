#!/usr/bin/env python3
"""
lrc14_far_order_leading_term_macmini_0620s5.py  (mac-mini-2026-06-20-S5)

REFRAMING E: the FAR-ORDER leading-term decomposition.

By THM-551, for a bounded core B the Newton far-expansion
    p0(B u F) = sum_{S subset F} Delta_S(B),   Delta_S(B) = sum_{T subset S}(-1)^{|S|-|T|} p0(B u T)
BEGINS at far-order s0 = max(0, 7-|B|) and (conjecturally) each successive order falls by ~1/7.

Decorrelated limit (THM-548 / HYP-2680):
    Phi_s(B) = 7^{-s} sum_{t=1}^{s} (-1)^{s-t} t! Stirling2(s,t) p_t(B),  p_t = meas{B misses exactly t sectors}.

This script, for a fixed true-wide row E = B u F (B bounded core, F far runners), computes:
  (A) The EXACT order-by-order sums
        D_s = sum_{|S|=s, S subset F} Delta_S(B)       (the actual far packets at order s)
        Ph_s = C(|F|,s)-weighted? NO -- Phi_s is per-subset; we sum it over the C(r,s) subsets of size s.
        R_s = D_s - C(r,s)*Phi_s(B)                    (the per-order RESIDUAL, reframing E's target)
      and the partial reconstruction  p0 ?= sum_s D_s  (exact identity check).
  (B) The geometric-decay test: ratio |D_s|/|D_{s-1}| and |R_s|/|R_{s-1}| vs 1/7.
  (C) The leading-order bound: is p0 <= (leading D_{s0}) + (geometric tail)? Report
        LEAD = D_{s0},  TAIL = sum_{s>s0} D_s,  and whether |TAIL| <= |LEAD|/6 (1/7-geom tail).
  (D) The "two sevens" check: at the dangerous rows |B| chosen so s0 = 7-|B| = 3, the leading
      order is THREE-far -> tie to 3+3+1 cube-root packet grading.

Uses the EXACT engine (rational measS7); F chosen DISSOCIATED (Sidon-like) so the decorrelated
limit is approached, and also a RESONANT F to see the residual blow up.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

# ---------------- exact engine ----------------
def sector_of(p): return int((p % 1) * 7)

def measS7(E):
    E = sorted(set(E)); bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7*e + 1): bps.add(F(m, 7*e))
    bps = sorted(bps); tot = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E)) == 7: tot += x1 - x0
    return tot

# miss-profile p_t(B) = meas{B misses exactly t sectors}
def miss_profile(B):
    B = sorted(set(B)); bps = set([F(0), F(1)])
    for e in B:
        if e == 0: continue
        for m in range(0, 7*e + 1): bps.add(F(m, 7*e))
    bps = sorted(bps); prof = {}
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        nmiss = 7 - len(set(sector_of(e*((x0+x1)/2)) for e in B))
        prof[nmiss] = prof.get(nmiss, F(0)) + (x1 - x0)
    return prof  # prof[t] = p_t(B)

# Stirling numbers of the second kind
def stirling2(n, k):
    if k == 0: return 1 if n == 0 else 0
    if k > n: return 0
    S = [[0]*(k+1) for _ in range(n+1)]
    S[0][0] = 1
    for i in range(1, n+1):
        for j in range(1, k+1):
            S[i][j] = j*S[i-1][j] + S[i-1][j-1]
    return S[n][k]

def fact(n):
    r = 1
    for i in range(2, n+1): r *= i
    return r

def Phi_s(B, s):
    """Phi_s(B) = 7^{-s} sum_{t=1}^{s} (-1)^{s-t} t! S(s,t) p_t(B)."""
    if s == 0: return measS7(B)  # Phi_0 = p0(B)
    prof = miss_profile(B)
    tot = F(0)
    for t in range(1, s+1):
        pt = prof.get(t, F(0))
        tot += F((-1)**(s-t)) * fact(t) * stirling2(s, t) * pt
    return tot / F(7)**s

# Delta_S(B) = sum_{T subset S} (-1)^{|S|-|T|} p0(B u T)
def Delta_S(B, S):
    S = list(S); tot = F(0); s = len(S)
    for r in range(s+1):
        for T in itertools.combinations(S, r):
            tot += F((-1)**(s-r)) * measS7(list(B)+list(T))
    return tot

def C(n, k):
    if k < 0 or k > n: return 0
    num = 1; den = 1
    for i in range(k): num *= (n-i); den *= (i+1)
    return num // den

# ---------------- main per-order analysis ----------------
def analyze(B, Fset, label):
    B = tuple(B); Fset = tuple(Fset); r = len(Fset)
    s0 = max(0, 7 - len(B))
    print("="*78)
    print(f"ROW: {label}   B={B} (|B|={len(B)})   F={Fset} (r={r})   leading far-order s0={s0}")
    p0_full = measS7(list(B)+list(Fset))
    print(f"  p0(B u F) [exact] = {p0_full} = {float(p0_full):.6f}")
    # per-order packets
    Dsum = {}  # D_s
    Rsum = {}  # R_s = D_s - C(r,s)*Phi_s
    for s in range(0, r+1):
        Ds = F(0)
        for S in itertools.combinations(Fset, s):
            Ds += Delta_S(B, S)
        Dsum[s] = Ds
        Phis = Phi_s(B, s)
        Rsum[s] = Ds - F(C(r, s)) * Phis
    # identity check: sum_s D_s == p0_full
    recon = sum(Dsum[s] for s in range(0, r+1))
    print(f"  IDENTITY  sum_s D_s = {recon}  (== p0?  {recon == p0_full})")
    print(f"  {'s':>2}{'D_s (order sum)':>22}{'C(r,s)*Phi_s':>22}{'R_s residual':>22}{'|D_s|/|D_{s-1}|':>16}")
    prevD = None
    for s in range(0, r+1):
        Ds = Dsum[s]; Rs = Rsum[s]; Phis = Phi_s(B, s); cPhi = F(C(r,s))*Phis
        ratio = (abs(float(Ds))/abs(float(prevD))) if (prevD not in (None, 0) and prevD != F(0)) else float('nan')
        print(f"  {s:>2}{float(Ds):>22.8f}{float(cPhi):>22.8f}{float(Rs):>22.8f}{ratio:>16.4f}")
        prevD = Ds
    # leading-order + geometric tail
    lead = Dsum[s0]
    tail = sum(Dsum[s] for s in range(s0+1, r+1))
    print(f"  LEAD (order s0={s0}) D_{s0} = {float(lead):.6f}")
    print(f"  TAIL (orders > s0)        = {float(tail):.6f}")
    if lead != 0:
        print(f"  |TAIL|/|LEAD| = {abs(float(tail))/abs(float(lead)):.4f}   (1/7-geom predicts <= 1/6 = 0.1667)")
    return Dsum, Rsum, p0_full

print("\n##### REFRAMING E: FAR-ORDER LEADING-TERM DECOMPOSITION #####\n")

# --- dangerous-row models. Goal: |B| with s0 = 7-|B| = 3  =>  |B| = 4 (the 3-far leading row). ---
# Also |B|=7 (s0=0, p0=p0(B)+corrections) and |B|=5 (s0=2).

# (1) |B|=4 core, r=3 DISSOCIATED far set (Sidon-like, large, coprime gaps)
analyze([0,1,2,3], [50,71,103], "|B|=4 s0=3, DISSOCIATED far (3-far leading)")

# (2) |B|=4 core, r=3 RESONANT far set (consecutive -> strong (1,-1) resonance)
analyze([0,1,2,3], [40,41,42], "|B|=4 s0=3, RESONANT far (consec)")

# (3) |B|=5 core, r=3, s0=2 (two-far leading)
analyze([0,1,2,4,8], [50,71,103], "|B|=5 s0=2, DISSOCIATED far (2-far leading)")

# (4) |B|=7 core (consec_7), r=3, s0=0  (p0(B)+corrections)
analyze([0,1,2,3,4,5,6], [50,71,103], "|B|=7 s0=0, DISSOCIATED far")

# (5) |B|=4 with r=4 to see one more order of decay (dissociated)
analyze([0,1,2,3], [50,71,103,149], "|B|=4 s0=3, r=4 DISSOCIATED (extra order)")

print("\nDONE.")
