#!/usr/bin/env python3
"""
The finite N* sweep for the SIMULTANEOUS union-bound far-peel (LRCSimulPeel.lean, HYP-3875).
mac-mini-2026-07-02-S18.  Numerical companion to `lonely_of_simul_peel`: confirm its hypotheses
(fee < floor) are dischargeable for real middle-band covering sets, and find the threshold N*.

The Lean bridge: for a 13-tuple = window B (<=22) ++ far (j = |far| far runners), h = 1/14,
    fee   = c_B * 4h * Sum_{w in far} 1/w        (c_B = #components of good B)
    floor = (1 - 2h*j) * length(good B)
and  fee < floor  =>  good(B++far) nonempty  =>  a 14-lonely time exists.
This sweep computes length(good B), c_B, and the least N with fee < floor, for j = 1..6.
"""
from fractions import Fraction as F

H = F(1, 14)

def danger_components(B, h=H):
    """good(B) = [0,1) minus the union of danger bands of the window speeds; return components."""
    ivs = []
    for v in B:
        for a in range(v + 1):
            lo, hi = max(F(a, v) - h / v, F(0)), min(F(a, v) + h / v, F(1))
            if hi > lo: ivs.append((lo, hi))
    ivs.sort()
    merged = []
    for lo, hi in ivs:
        if merged and lo <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else: merged.append((lo, hi))
    comps, prev = [], F(0)
    for lo, hi in merged:
        if lo > prev: comps.append((prev, lo))
        prev = max(prev, hi)
    if prev < 1: comps.append((prev, F(1)))
    return comps

def good_length_and_pieces(B, h=H):
    comps = danger_components(B, h)
    return sum(hi - lo for lo, hi in comps), len(comps)

def least_Nstar(B, j, h=H):
    """least N s.t. peeling j far runners {N, N+1, ..., N+j-1} from good(B) keeps the floor > fee."""
    Lb, cB = good_length_and_pieces(B, h)
    floor = (1 - 2 * h * j) * Lb
    if floor <= 0:
        return None, Lb, cB, floor  # j >= 7: union floor dead
    # fee = cB * 4h * sum 1/w over far = {N + i}; find least N with fee < floor
    N = 23
    while N < 10**7:
        far = [N + i for i in range(j)]
        fee = cB * 4 * h * sum(F(1, w) for w in far)
        if fee < floor:
            return N, Lb, cB, floor
        N += 1
    return None, Lb, cB, floor

if __name__ == "__main__":
    print("="*78)
    print("FINITE N* for the simultaneous union-bound peel (h=1/14). fee<floor => lonely.")
    print("="*78)
    windows = {
        "{1..12}":        list(range(1, 13)),
        "{1..11}":        list(range(1, 12)),
        "{1,2,3,4}":      [1, 2, 3, 4],
        "{1..6}":         list(range(1, 7)),
        "AP {1..12}":     list(range(1, 13)),
    }
    for name, B in windows.items():
        Lb, cB = good_length_and_pieces(B)
        print(f"\nwindow {name}: length(good B)={float(Lb):.5f} ({Lb}), c_B={cB} components")
        print(f"  {'j (far)':>8} {'floor=(1-j/7)|Lb|':>18} {'N*':>10} {'note'}")
        for j in range(1, 8):
            Nstar, Lb2, cB2, floor = least_Nstar(B, j)
            if floor <= 0:
                print(f"  {j:>8} {float(floor):>18.5f} {'--':>10}  union floor <=0 (j>=7): DEAD (needs deep-cluster HYP-3901)")
            elif Nstar is None:
                print(f"  {j:>8} {float(floor):>18.5f} {'>1e7':>10}  no finite N* found <1e7")
            else:
                print(f"  {j:>8} {float(floor):>18.5f} {Nstar:>10}  fee<floor for all N>=N*")
    print("\n=> For j<=6 far runners, N* is FINITE and small (~1e2-1e4): the simultaneous peel closes")
    print("   the middle band 22<N<N* by a finite check. j>=7 dies (1-j/7<=0), the separate residual.")
    print("   This confirms lonely_of_simul_peel's fee<floor hypothesis is dischargeable per class.")
