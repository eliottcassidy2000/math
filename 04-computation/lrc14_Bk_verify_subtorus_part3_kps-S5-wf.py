#!/usr/bin/env python3
"""
PART 3 adversarial verification: the DECISIVE THRESHOLD SPLIT and the
LOAD-BEARING claim '1/7 -> consecutive is the bounded-shape minimizer'.

kind-pasteur-2026-06-18-S5-wf

The angle's actual deliverable (since it concedes 2/7 B(k) is NOT closed and the
2/7 object is REFUTED) is:
  CLAIM-1/7:  For threshold 1/7, over bounded primitive shapes, the CONSECUTIVE AP
              {0,1,...,k-1} is the minimizer (no perforation helps), at k=7..10.
  CLAIM-2/7:  For threshold 2/7, perforation BEATS consecutive (mu_min < a(k)).

If CLAIM-1/7 holds, then (combined with scale-invariance: single run = consec)
it grounds the working 1/7 route's 'consecutive minimizes' structural fact.

We test BOTH exhaustively over bounded shapes (cap on spread), EXACT rationals.
A single perforated shape with mu_{1/7} < a17(k) REFUTES CLAIM-1/7.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import importlib.util, os

spec = importlib.util.spec_from_file_location(
    "vm", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "lrc14_Bk_verify_subtorusrelation_kps-S5-wf.py"))
vm = importlib.util.module_from_spec(spec); spec.loader.exec_module(vm)
mu_rig = vm.mu_rig


def prim(E):
    g = 0
    for e in E: g = gcd(g, e)
    return g == 1


def bounded_min(k, cap, thr):
    """exact min mu over primitive shapes {0}∪combo, combo subset of [1..cap], |E|=k."""
    best, bestE = None, None
    for combo in combinations(range(1, cap+1), k-1):
        E = (0,)+combo
        if not prim(E):
            continue
        v = mu_rig(list(E), thr)
        if best is None or v < best:
            best, bestE = v, E
    return best, bestE


if __name__ == "__main__":
    # canon 1/7 consecutive values (claimed consecutive = minimizer)
    a17 = {8:F(691,735), 9:F(247,294), 10:F(38,49), 11:F(1381,2205),
           12:F(13823,24255), 13:F(477,1078)}

    print("="*90)
    print("(A) recompute 1/7 consecutive a17(k) and confirm canon")
    print("="*90)
    for k in range(8, 14):
        v = mu_rig(list(range(k)), F(1,7))
        ok = (v == a17[k])
        print(f"  a17({k:2d}) = {str(v):>16s}  canon {str(a17[k]):>16s}  "
              f"{'MATCH' if ok else '*** MISMATCH ***'}")

    print("\n" + "="*90)
    print("(B) CLAIM-1/7: is consecutive the BOUNDED-SHAPE MINIMIZER for 1/7?")
    print("    exhaustive over spread cap; report min and whether it equals consec")
    print("="*90)
    # caps chosen so exhaustion is feasible; spread a few beyond k.
    plan = {7:(F(1,7),12), 8:(F(1,7),12), 9:(F(1,7),13), 10:(F(1,7),14)}
    for k,(thr,cap) in plan.items():
        b, E = bounded_min(k, cap, thr)
        consec = mu_rig(list(range(k)), thr)
        verdict = "CONSEC IS MIN" if b == consec else f"*** PERFORATION BEATS CONSEC: {E} ***"
        print(f"  k={k} cap={cap} thr=1/7: min={str(b):>14s} consec={str(consec):>14s}  {verdict}")

    print("\n" + "="*90)
    print("(C) CLAIM-2/7: perforation BEATS consecutive for 2/7 (mu_min<a(k))")
    print("    exhaustive over spread cap; report min")
    print("="*90)
    a27 = {7:F(83,210), 8:F(163,490), 9:F(277,980)}
    plan2 = {7:12, 8:12, 9:13}
    for k,cap in plan2.items():
        b, E = bounded_min(k, cap, F(2,7))
        consec = a27[k]
        verdict = ("PERFORATION BEATS CONSEC" if b < consec
                   else "consec is min (claim would be wrong)")
        print(f"  k={k} cap={cap} thr=2/7: min={str(b):>14s}={float(b):.5f} "
              f"consec={str(consec):>10s}={float(consec):.5f}  {verdict}  at {E}")
