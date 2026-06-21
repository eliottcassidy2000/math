#!/usr/bin/env python3
"""
lrc14_threadC_truemax_gate_macmini_0620.py   (mac-mini-2026-06-20, THREAD C)

THE CORRECTED S3 GATE for LRC(14).

PRIOR CLAIM (CLAUDE.md / HYP-2697): "consec_k MAXIMIZES measS7 over the finite family"
   -- VERIFIED only k=8,9,10. This script EXTENDS to k=11,12,13 and finds:

   * consec_k IS the maximizer for k = 8, 9, 10, 11.
   * consec_k is NOT the maximizer for k = 12, 13.
       k=12 maximizer: {0,1,...,10,12}     measS7 = 11381/17640 ~ 0.645182  (> consec 0.624114)
       k=13 maximizer: {0,1,...,10,12,13}  measS7 = 19492/28665 ~ 0.679993  (> consec 0.675928)
     The maximizer always keeps the consecutive prefix {0,...,10} and pushes the remaining
     offsets to the top of the window (skipping 11).

   * BUT the S3 GATE STILL HOLDS: the TRUE maximum of measS7(E) over the full span<=13
     family is < cap_k for EVERY k=8..13, with slack that GROWS in k:
        k :  true_max      cap_k       slack
        8 : 0.327211     0.381463    +0.054252
        9 : 0.416156     0.494256    +0.078099
       10 : 0.504478     0.604396    +0.099917
       11 : 0.581463     0.725275    +0.143812
       12 : 0.645182     0.857143    +0.211962
       13 : 0.679993     1.000000    +0.320007

   So the LRC(14) S3 closure is SAFE, but its JUSTIFICATION must be "true-max < cap" (a
   4095-set finite check), NOT "consec is the extremal shape" (which is FALSE for k>=12).

This is a full, exhaustive, exact-rational verification of the gate over the LRC(14) family.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)


def measS7_E(E):
    """measS7(E) exact via theta=7x breakpoints. E = offset set, 0 in E."""
    E = sorted(set(E))
    Enz = [e for e in E if e != 0]
    bps = set([F(0), F(7)])
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, e))
    bps = sorted(x for x in bps if 0 <= x <= 7)
    total = F(0)
    for i in range(len(bps) - 1):
        t0, t1 = bps[i], bps[i + 1]
        if t1 <= t0:
            continue
        tm = (t0 + t1) / 2
        if len(set(int(e * tm) % 7 for e in E)) == 7:
            total += t1 - t0
    return total / 7


def measGP(P):
    """cap building block: meas{x: all speeds p in P keep frac(px) >= 1/14 from 0 (and 1)}."""
    P = sorted(set(P))
    if not P:
        return F(1)
    bps = set([F(0), F(1)])
    th = F(1, 14)
    for p in P:
        for j in range(0, p + 1):
            for off in (F(1, 14), F(-1, 14)):
                v = (F(j) + off) / p
                if 0 <= v <= 1:
                    bps.add(v)
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for t in range(len(bps) - 1):
        x0, x1 = bps[t], bps[t + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        if all(min((p * xm) % 1, 1 - ((p * xm) % 1)) >= th for p in P):
            tot += x1 - x0
    return tot


CAP_P = {8: [1, 5, 7, 8, 9], 9: [1, 11, 12, 13], 10: [1, 12, 13],
         11: [1, 13], 12: [1], 13: []}


def main():
    print("=" * 96)
    print("CORRECTED S3 GATE: true max_{span<=13} measS7(E)  vs  cap_k,  exact, exhaustive")
    print("=" * 96)
    print(f"  {'k':>3} {'#sets':>6} {'consec=max?':>11} {'true_max':>16} {'maximizer':>32} "
          f"{'cap_k':>16} {'max<cap?':>9} {'slack':>9}")
    all_ok = True
    for k in range(8, 14):
        consec = tuple(range(k))
        cap = measGP(CAP_P[k])
        best = F(-1)
        bestE = None
        nsets = 0
        for combo in itertools.combinations(range(1, 14), k - 1):
            E = (0,) + combo
            nsets += 1
            v = measS7_E(E)
            if v > best:
                best = v
                bestE = E
        is_consec = (bestE == consec)
        ok = best < cap
        all_ok &= ok
        print(f"  {k:>3} {nsets:>6} {str(is_consec):>11} {str(best):>16} {str(bestE):>32} "
              f"{str(cap):>16} {str(ok):>9} {float(cap - best):>+9.5f}")
    print()
    print(f"  S3 GATE (true-max < cap_k) holds for ALL k=8..13 : {all_ok}")
    print("  consec is the maximizer ONLY for k=8,9,10,11; FALSE for k=12,13.")
    print("  Extremal shape (k>=12) = {0,...,10} + top-of-window fill (skips 11).")


if __name__ == "__main__":
    main()
