#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeB_convex_order_NE_opus_0620s1.py   (opus-2026-06-20-S1)

ROUTE B for HYP-2693/2607/2694 (the COMPRESSION/EXTREMALITY crux):
   prove  U4(E) = p_0 + p_5 + 5 p_6  <=  U4(consec_k)  for every bounded cluster E.
where N_E(x) = # of inner sectors {1..6} MISSED by orbit {frac(e_i x)}, p_t=meas{N=t},
and U4 = 1 - S_1 + S_2 - S_3 + S_4  (THM-556, PROVED identity).

Route B vehicle: a CONVEX / MAJORIZATION order on the distribution of N_E.

KEY ALGEBRAIC FACT (this script proves it numerically then we argue it):
  U4 weights phi = (phi_0..phi_6) = (1,0,0,0,0,1,5).  Second differences of phi:
     1, 0, 0, 1, 3   (all >= 0)  =>  phi is CONVEX on {0,..,6}.
  Hence  N_E <=_cx N_consec  =>  U4(E) <= U4(consec).
  Convex order requires EQUAL MEAN.  We test if E[N_E] is constant across E (it is NOT
  in general), so the real order is the INCREASING-CONVEX or a TAIL/cut-based order.

This script:
  (A) print exact p-vectors, E[N], Var[N], U4 for consec vs a large bank of bounded E.
  (B) test the three candidate stochastic orders against U4 ranking:
        - convex order   (equal mean required) via the "integrated survival twice" test
        - increasing-convex order   icx  (cut on E[(N-a)_+] for all a)
        - the EXACT cut functions that U4 is a positive combination of.
  (C) DECOMPOSE phi into a positive basis of "convex corner" functions and test
      whether consec dominates E on EACH basis cut function (a genuinely JOINT but
      cut-by-cut certificate).  If yes on all cuts -> majorization proof.
  (D) the difference-multiset majorization probe: does {e_i-e_j} of consec majorize
      that of any bounded E, and does U4 respect it?
"""
import sys, itertools
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact p-vector engine (breakpoints x=a/(7e)) ----------
def dist_p(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(F(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = (e*mid) % 1
            hit.add((v.numerator*7)//v.denominator)
        t = sum(1 for j in range(1,7) if j not in hit)
        p[t] += (hi-lo)
    return p

PHI = [F(1),F(0),F(0),F(0),F(0),F(1),F(5)]   # U4 weights
def U4(p): return sum(p[t]*PHI[t] for t in range(7))
def EN(p): return sum(t*p[t] for t in range(7))
def VarN(p):
    m = EN(p); return sum(t*t*p[t] for t in range(7)) - m*m

def second_diff(phi):
    return [phi[t-1]-2*phi[t]+phi[t+1] for t in range(1,6)]

# ---------- cut / survival functions ----------
def survival(p):  # G_a = P(N >= a), a=0..7 (G_0=1, G_7=0)
    return [sum(p[t] for t in range(a,7)) for a in range(8)]
def cut_above(p, a):  # E[(N-a)_+] = sum_{t>a}(t-a)p_t  -- icx test functions
    return sum((t-a)*p[t] for t in range(a+1,7))
def integrated_survival(p):  # H_a = sum_{b>=a} G_b = E[(N-a+1)_+]-ish; we use cut_above
    return [cut_above(p,a) for a in range(7)]

def consec(k): return list(range(k))

if __name__ == "__main__":
    print("="*78)
    print("ROUTE B: convex/cut order on N_E for U4 = p_0+p_5+5p_6")
    print("="*78)
    print("\nU4 weights phi =", [str(x) for x in PHI])
    print("second differences of phi (t=1..5):", [str(x) for x in second_diff(PHI)],
          " => phi CONVEX" if all(x>=0 for x in second_diff(PHI)) else " => NOT convex")
    print("phi monotone increasing?", all(PHI[t+1]>=PHI[t] for t in range(6)),
          " (=> NOT plain icx since it dips at t=1)")

    for k in (8, 9):
        print("\n" + "#"*70)
        print(f"# k={k}")
        print("#"*70)
        C = consec(k); pc = dist_p(C); Uc = U4(pc)
        print(f"consec={C}")
        print(f"  p = {[str(x) for x in pc]}")
        print(f"  E[N]={float(EN(pc)):.5f} Var={float(VarN(pc)):.5f} U4={Uc}={float(Uc):.6f}")
        print(f"  survival G_a (a=0..7) = {[str(x) for x in survival(pc)]}")
        print(f"  cut E[(N-a)_+] (a=0..6) = {[str(x) for x in integrated_survival(pc)]}")

        # ---- build the bank of bounded E (0 in E, |E|=k, span bounded) ----
        # use spread up to k+4 for thoroughness at k=8, k+3 at k=9
        span = k+4 if k==8 else k+3
        bank = []
        for rest in itertools.combinations(range(1, span+1), k-1):
            E = [0]+list(rest)
            bank.append(E)
        print(f"\n  bank size (span<= {span}): {len(bank)}")

        # (A) does any E beat consec on U4?
        beats = []
        # (B) order tests: count how often E violates each candidate order vs consec
        viol_cx = 0   # convex order N_E <=cx consec needs equal mean; count mean mismatches
        viol_icx_cut = 0  # cut_above(E,a) <= cut_above(consec,a) for ALL a? (icx upper)
        # we want consec to be the MAX in icx upper sense: cut_above(consec,a) >= cut_above(E,a)
        same_mean = 0
        max_EN = EN(pc); min_EN = EN(pc)
        worst_cut_a = None
        cutc = [cut_above(pc,a) for a in range(7)]
        survc = survival(pc)
        for E in bank:
            p = dist_p(E); U = U4(p); en = EN(p)
            if U > Uc + F(1,10**15):
                beats.append((E, U, p))
            if en == EN(pc): same_mean += 1
            max_EN = max(max_EN, en); min_EN = min(min_EN, en)
            # icx-upper: consec should dominate cut_above for every a
            cutE = [cut_above(p,a) for a in range(7)]
            if any(cutE[a] > cutc[a] + F(1,10**15) for a in range(7)):
                viol_icx_cut += 1
        print(f"\n  (A) shapes beating consec on U4: {len(beats)}")
        for (E,U,p) in beats[:8]:
            print(f"      {E}: U4={float(U):.6f} (> {float(Uc):.6f})  p={[float(x) for x in p]}")
        print(f"\n  (B) E[N] range over bank: [{float(min_EN):.5f}, {float(max_EN):.5f}]"
              f"  (consec E[N]={float(EN(pc)):.5f})")
        print(f"      shapes with EXACT same mean as consec: {same_mean}/{len(bank)}"
              f"   => convex order (equal-mean) {'APPLICABLE to all' if same_mean==len(bank) else 'NOT universally applicable'}")
        print(f"      shapes violating icx-upper cut dominance by consec "
              f"(some E[(N-a)_+] > consec): {viol_icx_cut}/{len(bank)}")
        print(f"        => if 0, consec dominates ALL increasing-convex functionals,"
              f" but U4 is NOT increasing so this is only suggestive.")
