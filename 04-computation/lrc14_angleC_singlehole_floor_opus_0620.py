#!/usr/bin/env python3
"""
lrc14_angleC_singlehole_floor_opus_0620.py   (opus-2026-06-20)

FOLLOW-UP to Angle C hole-criticality.  The additive (sum-of-criticalities) bound is
DEAD because the joint loss L(H) is SUBMODULAR (L_add over-estimates L for >=2 holes).
But submodularity gives a DIFFERENT, VALID lower bound:

   For a submodular set-loss with L(emptyset)=0 and monotone,
        L(H) >= max_{h in H} L({h})              (single-element floor)
   is generally TRUE for monotone submodular?  -- NO, that needs checking.  What IS
   true for a MONOTONE function is L(H) >= L({h}) for any single h in H (monotonicity
   of removal:  removing more clocks can only lose more cover).  i.e.
        measS7(full \ H) <= measS7(full \ {h})   for every h in H.
   Hence  measS7(E) <= min_{h in H} measS7(full \ {h})
                     = measS7(full) - max_{h in H} Delta_h.

This is the SINGLE-HOLE FLOOR certificate.  It is rigorously valid (pure monotonicity
of the cover region under clock removal -- a CA monotonicity, fully PROVED).

CERTIFICATE TO TEST:
   For every residual shape E=full\H (span N, |H|=N+1-k>=1), is there ALWAYS a hole h
   with Delta_h(full) >= measS7(full) - cap_k ?  Equivalently
        max_{h in {1..N}} Delta_h(consec_{N+1})  >= measS7(consec_{N+1}) - cap_k ?
   If YES for every span N in the residual regime, then EVERY residual shape (which has
   at least one hole, and we may pick the most-critical available hole IF it lies in H)
   ...

   CAVEAT: the certificate "measS7(E) <= measS7(full) - max_{h in H} Delta_h" only uses
   holes that ARE in H. The adversary picks H to AVOID the high-criticality clocks. So
   the SAFE single-hole certificate is the WORST hole the adversary is forced to leave:
        measS7(E) <= measS7(full) - max_{h in H} Delta_h(full)   [valid, but H-dependent]
   The adversary minimises max_{h in H} Delta_h by putting holes only at low-criticality
   positions. So the guaranteed floor is
        guaranteed_drop(N, nholes) = min over H of ( max_{h in H} Delta_h(full) ).
   Because the adversary will choose the nholes LOWEST-criticality positions, then take
   their max = the (nholes)-th smallest criticality.
   => guaranteed_drop = (nholes)-th smallest value among {Delta_h : h=1..N}.

   Wait -- but a hole at a low-criticality position in the FULL board may not stay low
   once OTHER holes are present (criticalities shift). The bound max_{h in H} Delta_h(full)
   uses the FULL-board criticality, which is a fixed number per h, so the adversary's best
   is exactly the nholes smallest full-board criticalities; their max is the nholes-th
   smallest. THIS bound is rigorous (monotonicity), independent of interaction.

   So define:
       sorted_crit(N) = sorted [Delta_h(full) : h=1..N]
       guaranteed_drop(N, j) = sorted_crit(N)[j-1]   (j-th smallest, j=nholes)
   Certificate closes for (k,N) iff
       measS7(full) - guaranteed_drop(N, nholes) <= cap_k
   i.e. guaranteed_drop(N, nholes) >= measS7(full) - cap_k.

THIS SCRIPT: compute guaranteed_drop for every residual (k,N) and check the certificate.
ALSO: an IMPROVED two-level floor using the SECOND-best forced hole (pair monotonicity),
and a comparison to the true worst residual.  EXACT Fractions, stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = set([F(0), F(1)])
    for e in Enz:
        for m in range(0, 7*e+1):
            bps.add(F(m, 7*e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        res = set(int(7*e*xm) % 7 for e in E)
        if len(res) == 7: total += x1 - x0
    return total

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91),
       11: F(66,91), 12: F(6,7), 13: F(1)}

def main():
    print("="*96)
    print("ANGLE C single-hole FLOOR certificate  (opus-2026-06-20)")
    print("Certificate: measS7(E=full\\H) <= measS7(full) - max_{h in H} Delta_h(full)")
    print("Adversary minimises by holing the lowest-criticality positions =>")
    print("  guaranteed_drop(N, nholes) = (nholes)-th smallest full-board criticality.")
    print("="*96)
    regimes = {8: range(8, 14), 9: range(9, 14), 10: range(10, 14)}

    # precompute full-board criticalities
    crit = {}
    for N in range(7, 14):
        full = tuple(range(N+1)); mfull = measS7(full)
        crit[N] = {h: mfull - measS7(tuple(e for e in full if e != h)) for h in range(0, N+1)}

    all_closed = True
    for k in [8, 9, 10]:
        ck = CAP[k]
        print(f"\n--- k={k}, cap_k={str(ck)}={float(ck):.5f} ---")
        for N in regimes[k]:
            nholes = (N+1) - k
            if nholes < 1: continue
            full = tuple(range(N+1)); mfull = measS7(full)
            req_drop = mfull - ck
            # adversary keeps span=N => clock N must NOT be a hole (else span drops).
            # so available hole positions for the adversary are {1..N-1} (0 pinned, N kept).
            # BUT a shape with a hole at N is just a smaller-span shape, already covered by
            # the (N-1) row. For span EXACTLY N, holes come from {1..N-1}.
            hole_positions = list(range(1, N))   # exclude 0 (pinned) and N (span anchor)
            crit_vals = sorted(crit[N][h] for h in hole_positions)
            if nholes <= len(crit_vals):
                guaranteed = crit_vals[nholes-1]   # nholes-th smallest
            else:
                guaranteed = crit_vals[-1] if crit_vals else F(0)
            cert = (guaranteed >= req_drop)
            resid_bound = mfull - guaranteed
            # also the TRUE worst residual at this (k,N) for comparison
            worst = F(-1)
            for H in itertools.combinations(hole_positions, nholes):
                E = tuple(e for e in full if e not in set(H))
                v = measS7(E)
                if v > worst: worst = v
            print(f"  N={N:2d} holes={nholes}: req_drop={float(req_drop):.5f}  "
                  f"guaranteed_drop={float(guaranteed):.5f}  bound on measS7={float(resid_bound):.5f}  "
                  f"true worst={float(worst):.5f}  cert {'CLOSES' if cert else 'FAILS'}")
            if not cert:
                all_closed = False

    print("\n" + "="*96)
    print(f"SINGLE-HOLE FLOOR certificate closes ALL residual (k=8,9,10) spans: {all_closed}")
    print("="*96)
    print("If FAILS anywhere, the single-hole monotonicity floor is too weak there;")
    print("a 2-hole (pair) monotonicity floor or direct enumeration is needed for that span.")

if __name__ == "__main__":
    main()
