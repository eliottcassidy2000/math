#!/usr/bin/env python3
"""
lrc14_far_order_residual_decay_macmini_0620s5.py  (mac-mini-2026-06-20-S5)

SHARPENING reframing E. The first pass showed: for the genuinely-dangerous rows
(|B| ~ 7, s0=0) the ORDER SUMS D_s do NOT decay geometrically -- D_1 ~ D_0. So
"leading order + 1/7 tail" on D_s is FALSE at the dangerous row.

The CORRECT object is the per-order RESIDUAL R_s = D_s - C(r,s)*Phi_s, because the
decorrelated main term P_r(B) = sum_s C(r,s)*Phi_s already accounts for the bulk
(THM-548). p0 = P_r(B) + sum_s R_s. So reframing E is really:
    p0 = P_r(B)  +  sum_{s>=1} R_s,
and we test whether the RESIDUALS R_s decay geometrically (~1/7 per order), which
would let us bound the correction sum_s R_s by its leading residual + 1/7 tail.

This script, for the dangerous |B|=7 row with a DISSOCIATED far set at increasing
scale, computes R_s and P_r(B), tracks |R_s|/|R_{s-1}|, and checks
    |sum_{s>=1} R_s| <= |R_1| * (1/(1-1/7)) = |R_1|*7/6  ?
We ALSO scan multiple dissociated F at different scales to see if R_s -> 0 as F->inf
(the Fatou limit: residuals should VANISH for dissociated far runners).
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

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
    return prof
def stirling2(n, k):
    if k == 0: return 1 if n == 0 else 0
    if k > n: return 0
    S = [[0]*(k+1) for _ in range(n+1)]; S[0][0] = 1
    for i in range(1, n+1):
        for j in range(1, k+1): S[i][j] = j*S[i-1][j] + S[i-1][j-1]
    return S[n][k]
def fact(n):
    r = 1
    for i in range(2, n+1): r *= i
    return r
def Phi_s(B, s):
    if s == 0: return measS7(B)
    prof = miss_profile(B); tot = F(0)
    for t in range(1, s+1):
        tot += F((-1)**(s-t)) * fact(t) * stirling2(s, t) * prof.get(t, F(0))
    return tot / F(7)**s
def Delta_S(B, S):
    S = list(S); tot = F(0); s = len(S)
    for r in range(s+1):
        for T in itertools.combinations(S, r):
            tot += F((-1)**(s-r)) * measS7(list(B)+list(T))
    return tot
def Cnk(n, k):
    if k < 0 or k > n: return 0
    num=1; den=1
    for i in range(k): num*=(n-i); den*=(i+1)
    return num//den

def P_r(B, r):
    """Fully-decorrelated main term P_r(B) = sum_{t=0}^6 prof_t(B) c_t(r),
       c_t(r) = sum_i (-1)^i C(t,i)(1-i/7)^r."""
    prof = miss_profile(B); tot = F(0)
    for t in range(0, 7):
        pt = prof.get(t, F(0))
        ct = F(0)
        for i in range(0, t+1):
            ct += F((-1)**i) * Cnk(t, i) * F(7-i, 7)**r
        tot += pt * ct
    return tot

def analyze_residuals(B, Fset, label):
    B = tuple(B); Fset = tuple(Fset); r = len(Fset)
    print("-"*78)
    print(f"{label}: B={B} (|B|={len(B)})  F={Fset}")
    p0_full = measS7(list(B)+list(Fset))
    Pr = P_r(B, r)
    print(f"  p0(BuF)={float(p0_full):.6f}   P_r(B)[decorr main]={float(Pr):.6f}   corr=p0-P_r={float(p0_full-Pr):+.6f}")
    Rsum = {}
    for s in range(0, r+1):
        Ds = F(0)
        for S in itertools.combinations(Fset, s):
            Ds += Delta_S(B, S)
        Rsum[s] = Ds - F(Cnk(r, s)) * Phi_s(B, s)
    # consistency: sum_s R_s == p0 - P_r  (since sum_s C(r,s)Phi_s == P_r)
    sumR = sum(Rsum[s] for s in range(0, r+1))
    print(f"  sum_s R_s = {float(sumR):+.6f}   (should equal p0-P_r = {float(p0_full-Pr):+.6f}; match={sumR==p0_full-Pr})")
    print(f"  {'s':>2}{'R_s':>16}{'|R_s|/|R_{s-1}|':>18}")
    prev = None
    for s in range(1, r+1):
        Rs = Rsum[s]
        ratio = (abs(float(Rs))/abs(float(prev))) if (prev not in (None,) and prev != F(0)) else float('nan')
        print(f"  {s:>2}{float(Rs):>16.8f}{ratio:>18.4f}")
        prev = Rs
    R1 = Rsum.get(1, F(0))
    if R1 != 0:
        bound = abs(float(R1))*7/6
        print(f"  |sum_{{s>=1}} R_s| = {abs(float(sumR)):.6f}  vs  |R_1|*7/6 (geom-tail bound) = {bound:.6f}  -> {'OK' if abs(float(sumR))<=bound+1e-12 else 'VIOLATED'}")
    return p0_full, Pr, sumR

print("\n##### REFRAMING E (sharpened): per-order RESIDUAL decay #####\n")

print("### Dangerous |B|=7 (consec_7), dissociated F at growing scale: does corr->0 (Fatou)? ###")
for scale in [1,2,4,8]:
    F3 = (50*scale+1, 71*scale, 103*scale-1)  # keep coprime-ish, grow scale
    analyze_residuals([0,1,2,3,4,5,6], F3, f"scale={scale}")

print("\n### |B|=6, dissociated, growing scale ###")
for scale in [1,2,4]:
    F3 = (50*scale+1, 71*scale, 103*scale-1)
    analyze_residuals([0,1,2,3,4,5], F3, f"|B|=6 scale={scale}")

print("\n### |B|=7 RESONANT far (consec) -- residuals should NOT vanish ###")
for u in [40, 80, 160]:
    analyze_residuals([0,1,2,3,4,5,6], (u,u+1,u+2), f"consec u={u}")

print("\n### |B|=4 dissociated, r=4, residual decay (the 3-far leading dangerous shape) ###")
for scale in [1,2,4]:
    F4 = (50*scale+1, 71*scale, 103*scale-1, 149*scale)
    analyze_residuals([0,1,2,3], F4, f"|B|=4 r=4 scale={scale}")

print("\nDONE.")
