#!/usr/bin/env python3
"""
lrc14_angleC_worstshape_structure_opus_0620.py   (opus-2026-06-20)

Angle C, structural closure attempt.  The compaction/additive routes are DEAD (measS7
non-monotone under local moves).  But the EXHAUSTIVE residual data shows a clean fact:

   For fixed k, as span N grows, the worst (max-measS7) k-subset's value DROPS
   monotonically below consec_k = the N=k-1 value.

If we can prove the SINGLE statement
   W_k(N) := max{ measS7(E) : |E|=k, 0 in E, span(E)=N }  satisfies  W_k(N) <= consec_k
for all N >= k-1 (with equality only at N=k-1), AND W_k(N) -> 0, then combined with the
wide-span smallness bound we close LRC(14) for k=8,9,10.

This script:
 (1) tabulates W_k(N) and its argmax shape for N=k-1 .. up to a large box, EXACT.
 (2) tests whether W_k(N) is itself MONOTONE DECREASING in N for N>=k-1 (a 'span peels off'
     statement). If monotone, then proving W_k(k-1)=consec_k dominates W_k(k) suffices to
     chain (but need monotone, which we test honestly).
 (3) identifies the argmax SHAPE FAMILY (is it always 'consec block + isolated tail'?),
     which would give a 1-parameter family to bound in closed form.
 (4) HONEST: reports where W_k(N) is NOT monotone (the C2 phenomenon at the worst-case level)
     and whether the global sup over N is still consec_k.

EXACT Fractions; stdlib only.
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

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

def describe_shape(E):
    """Return a compact gap-signature: list of gaps between consecutive kept clocks."""
    E = sorted(E)
    return tuple(E[i+1]-E[i] for i in range(len(E)-1))

def main():
    print("="*96)
    print("ANGLE C worst-shape structure W_k(N)=max measS7 over span-N k-subsets  (opus-0620)")
    print("="*96)
    for k in [8, 9, 10]:
        ck = CAP[k]
        consec = measS7(tuple(range(k)))
        print(f"\n--- k={k}  consec_k(=W_k({k-1}))={float(consec):.6f}  cap_k={float(ck):.5f} ---")
        Wprev = None
        Wmax_over_N = F(-1); arg = None
        Nbox = 16
        for N in range(k-1, Nbox+1):
            # all k-subsets of {0..N} containing 0 and N (span exactly N)
            best = F(-1); bestE = None
            inner = N-1  # positions strictly between 0 and N
            need = k-2   # we already fix 0 and N
            if need < 0 or need > inner:
                continue
            for body in itertools.combinations(range(1, N), need):
                E = (0,)+body+(N,)
                v = measS7(E)
                if v > best: best, bestE = v, E
            mono = "" if Wprev is None else ("  (DOWN)" if best < Wprev else ("  (UP!)" if best > Wprev else "  (=)"))
            sig = describe_shape(bestE)
            print(f"  N={N:2d}: W_k(N)={float(best):.6f}  argmax={bestE}  gaps={sig}{mono}")
            Wprev = best
            if best > Wmax_over_N: Wmax_over_N, arg = best, bestE
        print(f"  >>> sup_N W_k(N) over N in [{k-1},{Nbox}] = {float(Wmax_over_N):.6f} at {arg}; "
              f"consec is the sup: {Wmax_over_N==consec}  (<=cap: {Wmax_over_N<=ck})")

    print("\n" + "="*96)
    print("Honest read: see whether W_k(N) is monotone-decreasing (clean peel) or has UP! steps.")
    print("="*96)

if __name__ == "__main__":
    main()
