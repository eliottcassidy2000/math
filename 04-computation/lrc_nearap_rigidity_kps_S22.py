#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S22: THE NEAR-AP RIGIDITY MECHANISM (the isolated hard
kernel of (G)) through the roots-of-unity / discrepancy lens.

State (fleet): mac-mini HYP-4402 reduced (G) to SINGLE-SCALE near-AP covering
families; opus HYP-4013/4074 the DISCREPANCY INVERSION (AP = 13th roots of unity
at t=1/13 = min star-discrepancy = min M; Fibonacci = max, loose).  The open
pieces: strict lift-rigidity (M=1/13 => dilated AP) + the density floor.

THE ROOTS-OF-UNITY PICTURE: at t=1/13 the AP {1..12} sits at the nonzero 13th
roots of unity, margin 1/13, achieved by the TWO critical runners r=1 (at +1/13)
and r=12 (at -1/13).  As t moves off 1/13, runner r moves at rate v_r; the AP is
a BALANCED saddle (runner 1 wants t up, runner 12 wants t down) so M stays 1/13.
A LIFT (v_r = r + 13 c_r) changes the rates and BREAKS the balance => M > 1/13.

We test: (a) is the AP the UNIQUE M=1/13 family among all lifts (strict rigidity)?
(b) what is the infimum of M over nonzero lifts, and its mechanism? (c) does the
balance-breaking give a CLEAN sufficient condition (both critical runners' lifts
push M up)?  (d) the density-floor: min M over non-AP near-AP covering families.
"""
from fractions import Fraction
from itertools import product
import numpy as np

LO = Fraction(1, 13)   # AP value
HI = Fraction(2, 25)   # gap upper edge

def M_exact(v, Qcap=None):
    S = int(sum(abs(x) for x in v)); Q = Qcap if Qcap else 4 * S
    va = np.array(v, dtype=np.int64); bn, bd = 0, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q, dtype=np.int64); r = np.outer(va, a) % q
        d = np.minimum(r, q - r); bq = int(d.min(axis=0).max())
        if bq * bd > bn * q: bn, bd = bq, q
    return Fraction(bn, bd)

AP = list(range(1, 13))

print("=== (a) single-runner lifts of the AP: {1..12} with runner j -> j+13 ===", flush=True)
print("  (mac-mini: all have M >= 2/23; verifying + the mechanism)", flush=True)
minM_single = (Fraction(1), None)
for j in range(12):
    v = AP.copy(); v[j] = v[j] + 13
    M = M_exact(v)
    if M < minM_single[0]: minM_single = (M, v)
    print(f"  lift runner {j+1} (-> {v[j]}): M = {M} = {float(M):.5f}  {'STILL 1/13!' if M==LO else ('IN GAP' if LO<M<HI else 'above')}", flush=True)
print(f"  min single-lift M = {minM_single[0]} = {float(minM_single[0]):.5f} at {minM_single[1]}", flush=True)

print(flush=True)
print("=== (b) all lifts in {0,13} per runner (block structure), min M over NONZERO lifts ===", flush=True)
# each runner either stays r or lifts to r+13; 2^12 combos, find min M over nonzero
minM_block = (Fraction(1), None)
ap_val = M_exact(AP)
n_at_1_13 = 0; n_in_gap = 0; total = 0
for mask in range(1, 4096):   # nonzero lifts
    v = [AP[i] + (13 if (mask >> i) & 1 else 0) for i in range(12)]
    M = M_exact(v)
    total += 1
    if M == LO: n_at_1_13 += 1
    if LO < M < HI: n_in_gap += 1
    if M < minM_block[0]: minM_block = (M, v, mask)
print(f"  AP itself M = {ap_val}", flush=True)
print(f"  over {total} nonzero {{0,13}}-lifts: {n_at_1_13} still at 1/13, {n_in_gap} IN OPEN GAP", flush=True)
print(f"  min nonzero-lift M = {minM_block[0]} = {float(minM_block[0]):.5f}  (infimum should be 2/25 at block lift)", flush=True)
print(f"     minimizer bits={bin(minM_block[2])}, v={minM_block[1]}", flush=True)
if n_at_1_13 == 0:
    print("  => STRICT RIGIDITY holds over {0,13}-lifts: NO nonzero lift stays at 1/13.", flush=True)

print(flush=True)
print("=== (c) balance-breaking: the two critical runners (r=1, r=12) and lift direction ===", flush=True)
# lift only runner 1 and/or runner 12, see when M jumps
for (a1, a12) in [(1,0),(0,1),(1,1),(2,0),(0,2),(1,2),(2,1),(2,2)]:
    v = AP.copy(); v[0] += 13*a1; v[11] += 13*a12
    M = M_exact(v)
    print(f"  lift r=1 by {a1}*13, r=12 by {a12}*13: v[0]={v[0]}, v[11]={v[11]}, M = {float(M):.5f}", flush=True)

print(flush=True)
print("=== (d) density floor: min M over near-AP lifts (each runner in {r, r+13, r+26}) ===", flush=True)
# broader lift space, find the true infimum M over nonzero (non-dilation) lifts
minM_wide = (Fraction(1), None)
in_gap_wide = 0
import random
random.seed(2026070622)
for _ in range(30000):
    v = [AP[i] + 13*random.randint(0,3) for i in range(12)]
    if v == AP: continue
    # skip pure dilations d*{1..12}? check
    M = M_exact(v)
    if LO < M < HI:
        in_gap_wide += 1
        if in_gap_wide <= 5:
            print(f"  *** GAP MEMBER: v={v} M={float(M):.5f} ***", flush=True)
    if M < minM_wide[0] and M > LO:
        minM_wide = (M, v)
print(f"  over 30000 random near-AP lifts: {in_gap_wide} in the OPEN GAP", flush=True)
print(f"  min M (> 1/13) found = {minM_wide[0]} = {float(minM_wide[0]):.5f}", flush=True)
print(flush=True)
print("VERDICT: if 0 in-gap over all lift probes AND min-nonzero-lift = 2/25 => the density", flush=True)
print("floor is EXACTLY the 1/13->2/25 jump (strict rigidity + floor at the block lift).", flush=True)
