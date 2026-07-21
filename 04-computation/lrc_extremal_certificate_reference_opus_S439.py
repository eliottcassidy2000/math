#!/usr/bin/env python3
"""
REFERENCE COMPUTATION for the LRC-extremal native_decide Lean certificate (opus-S439).
The WOWII-103 template applied to LRC(14): certify M(V) = max_t min_i ||v_i t|| exactly, as a
FINITE decidable check over a critical grid of rationals -- the spec the Lean native_decide mirrors.

M(V) = max over t in [0,1) of min over v in V of dist(v*t, Z).  ||x|| = distance to nearest integer.
The piecewise-linear function t -> min_i ||v_i t|| has all its local maxima at "beat" times:
t = a/D with D in {v_i+v_j, |v_i-v_j|, 2*v_i}.  So M(V) = max over that FINITE grid -- exact, rational.
(This is the finite-attainment fact behind THM-401's pair-sum-time pinch, modulus 2n-1.)

We (1) confirm M({1..13}) = M({1..11,13,24}) = 1/14 on the beat grid, (2) cross-check against a fine
continuous scan, (3) emit the exact grid data + argmax for the Lean file, (4) report grid size
(native_decide feasibility).
"""
from fractions import Fraction as F
from itertools import combinations

def dist_to_int(x):
    """||x|| = distance to nearest integer, exact for Fraction x."""
    f = x - (x.numerator // x.denominator)   # frac part in [0,1)
    return min(f, 1 - f)

def gmin(V, t):
    return min(dist_to_int(v * t) for v in V)

def beat_denoms(V):
    D = set()
    for v in V:
        D.add(2 * v)
    for a, b in combinations(V, 2):
        D.add(a + b)
        D.add(abs(a - b))
    D.discard(0)
    return sorted(D)

def M_exact(V):
    """exact M(V) over the beat grid + all argmax t."""
    best = F(0); args = []
    for D in beat_denoms(V):
        for a in range(0, D):        # t = a/D in [0,1)
            t = F(a, D)
            g = gmin(V, t)
            if g > best:
                best = g; args = [t]
            elif g == best and t not in args:
                args.append(t)
    return best, sorted(set(args))

def M_scan(V, N=200000):
    """continuous cross-check: max over a fine uniform grid (float)."""
    best = 0.0
    Vf = [float(v) for v in V]
    for k in range(1, N):
        t = k / N
        g = min(abs(v*t - round(v*t)) for v in Vf)
        if g > best: best = g
    return best

FAMILIES = {
    "{1..13}  (the classical tight LRC(14) family)": list(range(1, 14)),
    "{1..11,13,24}  (leaf-inflation of {1..13}, 12->24)": list(range(1,12)) + [13, 24],
    "{1..11,13,36}  (the (1/14,3/41) gap witness)": list(range(1,12)) + [13, 36],
}

for name, V in FAMILIES.items():
    M, args = M_exact(V)
    grid = sum(D for D in beat_denoms(V))          # candidate-t count (<= this, before dedup)
    sc = M_scan(V)
    print(f"\n{name}")
    print(f"   V = {V}")
    print(f"   M(V) exact (beat grid) = {M}  = {float(M):.6f}")
    print(f"   continuous scan max    ~ {sc:.6f}   (cross-check)")
    print(f"   argmax t (first few)   = {args[:6]}{' ...' if len(args)>6 else ''}   ({len(args)} total)")
    print(f"   beat denominators      = {beat_denoms(V)}")
    print(f"   candidate-t upper bound (sum of denoms) = {grid}  -> native_decide-feasible")

# ---- emit Lean-ready data for the two 1/14 extremals ----
print("\n" + "="*66)
print("LEAN-READY: the two extremals both certify M = 1/14")
for name, V in list(FAMILIES.items())[:2]:
    M, args = M_exact(V)
    assert M == F(1,14), f"expected 1/14, got {M}"
    print(f"  V={V}:  M={M}  argmax includes t=1/14? {F(1,14) in args}")
print("\nBoth reduce to a finite max-of-min over rationals => native_decide certificate is a")
print("straight port: define V:List Nat, the beat grid, dist_to_int on Rat, and `by native_decide`.")
