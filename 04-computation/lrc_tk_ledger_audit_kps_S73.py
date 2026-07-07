#!/usr/bin/env python3
r"""
lrc_tk_ledger_audit_kps_S73.py   (kind-pasteur-2026-07-07-S73, HYP-5147 part d)

INDEPENDENT RECOMPUTATION of the (A') union-bound bars T_k, k=8..13.

Everything on the (A') critical path consumes the per-k thresholds
    T_k = m_P + 1 - min_{|P| = 13-k} meas(G_P),
quoted as {0.6185, 0.5057, 0.3956, 0.2747, 0.1429, 0.0565} (monad-explorer-S1, boxeph-S1).
This audit recomputes them EXACTLY (rational arithmetic, no reuse of monad's script logic
beyond the G_P definition), and records the argmin slow parts P.

G_P = {x in [0,1) : ||p x|| >= 1/14 for all p in P}   (THM-527/530 skeleton convention)
m_P = 14249/252252 (THM-530).
"""
from fractions import Fraction as F
from itertools import combinations

M_P = F(14249, 252252)

def gp_measure_exact(P):
    """meas(G_P) exactly: complement = union of intervals around j/p of radius 1/(14p)."""
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            lo, hi = F(j, p) - w, F(j, p) + w
            bad.append((max(lo, F(0)), min(hi, F(1))))
    bad.sort()
    tot = F(0); cur_lo = cur_hi = None
    for lo, hi in bad:
        if cur_lo is None:
            cur_lo, cur_hi = lo, hi
        elif lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            tot += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    if cur_lo is not None:
        tot += cur_hi - cur_lo
    return 1 - tot

print("=" * 88)
print("T_k LEDGER AUDIT: T_k = m_P + 1 - min_{|P|=13-k} meas(G_P), exact rationals")
print(f"m_P = {M_P} = {float(M_P):.6f}")
print("=" * 88)
QUOTED = {8: 0.6185, 9: 0.5057, 10: 0.3956, 11: 0.2747, 12: 0.1429, 13: 0.0565}
print(f"{'k':>3} {'|P|':>4} {'argmin P':>22} {'meas(G_P)':>24} {'T_k exact':>24} {'float':>8} {'quoted':>8} {'ok?':>4}")
for k in range(13, 7, -1):
    s = 13 - k
    if s == 0:
        best, bestP = F(1), ()
    else:
        best, bestP = F(2), None
        for P in combinations(range(1, 14), s):
            m = gp_measure_exact(P)
            if m < best:
                best, bestP = m, P
    Tk = M_P + 1 - best
    ok = abs(float(Tk) - QUOTED[k]) < 5e-4
    print(f"{k:>3} {s:>4} {str(bestP):>22} {str(best):>24} {str(Tk):>24} {float(Tk):8.4f} {QUOTED[k]:8.4f} {str(ok):>4}")

print()
print("FINDING: quoted fleet 'T_k' = 1 - meas(G_P) = the POSITIVITY bars.  The Lean skeleton's")
print("hlarge (LRCFourteenSkeleton.lean: witness_floor_from_cluster_cases) demands")
print("witnessMP <= witnessG2(shape) for EVERY 8<=k<=13 shape -- the QUANTITATIVE floor m_P.")
print("Via the union bound rho* >= meas(G_P) + mu - 1 >= m_P, the honest per-k tail bar is")
print("T_k = m_P + 1 - meas(G_P): exactly m_P = 0.0565 HIGHER at every k with |P|>=1.")
print("(MISTAKE-118 bar drift recurring inside the (A') ledger itself.)")
print()
print("=" * 88)
print("DOWNSTREAM IMPACT (fleet numbers vs honest bars)")
print("=" * 88)
onean  = {8: 0.602, 9: 0.511, 10: 0.434, 11: 0.368, 12: 0.321, 13: 0.281}   # boxeph-S1 1-anchor inf
twoan  = {8: 0.766, 9: 0.685, 10: 0.570, 11: 0.487, 12: 0.421, 13: 0.360}   # boxeph-S1 2-anchor inf
decorr = {8: 0.799, 10: 0.593, 13: 0.365}                                    # kps-S69 PA_2^inf limits
mins = {13: F(1), 12: F(6,7), 11: F(66,91), 10: F(55,91), 9: F(1979,4004), 8: F(2243,5880)}
print(f"{'k':>3} {'honest T_k':>11} {'1-anchor':>9} {'ok?':>4} {'2-anchor':>9} {'ok?':>4} {'margin':>7}")
for k in range(8, 14):
    Tk = float(M_P + 1 - mins[k])
    ok1 = onean[k] >= Tk; ok2 = twoan[k] >= Tk
    print(f"{k:>3} {Tk:>11.4f} {onean[k]:>9.3f} {str(ok1):>4} {twoan[k]:>9.3f} {str(ok2):>4} {twoan[k]-Tk:>7.3f}")
print()
print("=> 1-anchor route now FAILS k=8,9,10 (was: only k=8).  2-anchor still discharges all")
print("   k=8..13 but the k=8 margin shrinks 0.148 -> 0.091.  All (A')-side margin claims")
print("   (kps-S68 bites, S69/S70 margins) must be re-measured against the honest bars.")
