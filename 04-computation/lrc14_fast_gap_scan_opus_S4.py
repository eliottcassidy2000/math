#!/usr/bin/env python3
"""
lrc14_fast_gap_scan_opus_S4.py    opus-2026-07-23-S4
Fast EXACT gap(V)=max_tau min_v ||v tau|| for 13-speed configs, via integer arithmetic.
The max of the piecewise-linear g is at a breakpoint tau=k/d with
   d in D = {2v} u {|v_i-v_j|} u {v_i+v_j},
so gap = max_{d in D} max_{k=1..d-1} min_v min(vk mod d, d-vk mod d) / d.   EXACT (a rational).
"""
import numpy as np
from math import gcd
from functools import reduce
from fractions import Fraction as Fr

def gap_exact(V):
    """exact gap as Fraction, plus (best_d,best_k) so tau*=k/d."""
    n = len(V); Va = np.array(V, dtype=np.int64)
    D = set()
    for v in V: D.add(2 * v)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = V[i], V[j]
            if a != b: D.add(abs(a - b))
            D.add(a + b)
    best_num, best_den, bk = 0, 1, None
    for d in D:
        if d < 2: continue
        K = np.arange(1, d, dtype=np.int64)
        R = np.outer(Va, K) % d                 # (n, d-1)
        Dm = np.minimum(R, d - R)
        m = Dm.min(axis=0)                      # min over speeds, per k
        idx = int(m.argmax()); val = int(m[idx])
        if val * best_den > best_num * d:       # val/d > best_num/best_den
            best_num, best_den, bk = val, d, int(K[idx])
    return Fr(best_num, best_den), (best_den, bk)

def binding(V, tau):
    g = min(min((v * tau.numerator) % tau.denominator,
                tau.denominator - (v * tau.numerator) % tau.denominator) for v in V)
    return [v for v in V
            if min((v * tau.numerator) % tau.denominator,
                   tau.denominator - (v * tau.numerator) % tau.denominator) == g]

if __name__ == "__main__":
    FL = Fr(1, 14)
    print("VALIDATION against CONSTANTS-INDEX known 13-speed values:")
    tests = {
        "AP {1..13} (expect 1/14)":        list(range(1, 14)),
        "GW {1..11,13,24} (expect 1/14)":  list(range(1, 12)) + [13, 24],
        "{1..11,13,36} (expect 3/41)":     list(range(1, 12)) + [13, 36],
        "{1..12,26} (expect 2/27)":        list(range(1, 13)) + [26],
        "{1..12,39} (expect 3/40)":        list(range(1, 13)) + [39],
        "{1..11,13,48} (expect 4/53)":     list(range(1, 12)) + [13, 48],
        "{1..12,182} (expect 14/183)":     list(range(1, 13)) + [182],
        "{1..12,14} (expect 1/13)":        list(range(1, 13)) + [14],
    }
    ok = True
    for nm, V in tests.items():
        g, (d, k) = gap_exact(V)
        tau = Fr(k, d)
        b = binding(V, tau)
        print(f"  {nm:34s} -> gap={str(g):>8}  tau*={str(tau):>8}  binding={b}")
    print("\nAll should match the CONSTANTS-INDEX attained values.")
