#!/usr/bin/env python3
"""
parity_proof_s20cn.py -- kind-pasteur-2026-03-22-S20cn

THEOREM: Cross-orbit count F_n = T_n - SL_n is ALWAYS even.

Proof: The map Phi: (T, e) -> (T xor e, e) is an involution on
(tournament, arc) pairs that commutes with S_n. It maps cross-pairs
to cross-pairs (since T xor e in different class iff so). Since Phi
commutes with S_n, it induces an involution on orbits. A cross-orbit
O maps to Phi(O) != O (because no element of O satisfies T xor e = T'
with T' in same class as T). So cross-orbits come in pairs, and F is even.

CONSEQUENCE: T_n - SL_n is even. So SL_n and T_n have the same parity.

Verification: T_n is always even (from the Burnside formula, the identity
contributes 2^m * m which is always even).

So SL_n MUST be even. If our formula gives odd SL, it's WRONG.

Author: kind-pasteur-2026-03-22-S20cn
"""
import sys
from math import comb, factorial, gcd
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def partitions_odd(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0: yield []; return
    start = min(max_part, n)
    if start % 2 == 0: start -= 1
    for k in range(start, 0, -2):
        for rest in partitions_odd(n - k, k):
            yield [k] + rest

def ct_count(parts, n):
    mult = Counter(parts)
    d = 1
    for c, a in mult.items():
        d *= (c ** a) * factorial(a)
    return factorial(n) // d

def p_val(parts):
    r = len(parts)
    v = sum(gcd(parts[i], parts[j]) for i in range(r) for j in range(i+1, r))
    v += sum((c - 1) // 2 for c in parts)
    return v

print("=" * 70)
print("  PARITY ANALYSIS: T_n and SL_n must have same parity")
print("=" * 70)

# Known SL from computation
known_SL = {3: 2, 4: 6, 5: 16, 6: 58}

# NTE recurrence
NTE = {3: 0, 4: 48}
for nn in range(5, 20):
    NTE[nn] = comb(nn, 2) * NTE[nn-1]

print(f"\n  {'n':>3s} {'T_n':>12s} {'T_n%2':>5s} {'SL_twin':>8s} {'SL%2':>5s} {'Parity Match':>13s}")

for n in range(3, 16):
    m = comb(n, 2)
    nfact = factorial(n)

    # T_n
    T_sum = 0
    for parts in partitions_odd(n):
        ct = ct_count(parts, n)
        p = p_val(parts)
        f = parts.count(1)
        fix = 1 << p
        T_sum += ct * fix * comb(f, 2)
    T_n = T_sum // nfact

    # SL_n from twin formula
    twins = m * (1 << (m - n + 2))
    self_wt = twins + NTE[n]

    twin_corr = 0
    for parts in partitions_odd(n):
        if parts == [1] * n: continue
        ct = ct_count(parts, n)
        f = parts.count(1)
        if f < 2: continue
        p = p_val(parts)
        fix = 1 << p
        r = len(parts)
        twin_corr += ct * fix * comb(f, 2) / (1 << (r - 2))

    SL_numer = self_wt + twin_corr
    SL = SL_numer / nfact
    SL_int = int(round(SL))

    t_parity = T_n % 2
    sl_parity = SL_int % 2
    match = "OK" if t_parity == sl_parity else "*** FAIL ***"

    known = known_SL.get(n)
    known_str = f" (known={known})" if known else ""

    print(f"  {n:>3d} {T_n:>12d} {t_parity:>5d} {SL_int:>8d} {sl_parity:>5d} {match:>13s}{known_str}")

print(f"""
  THEOREM: F_n = T_n - SL_n is always even.
  Proof: The flip involution Phi(T,e) = (T+e, e) pairs cross-orbits.

  CONSEQUENCE: If our SL formula gives wrong parity,
  either the NTE recurrence or the twin formula is incomplete.

  At n=8: T=188288 (even), SL_formula=3145 (odd) -> PARITY VIOLATION!
  This means the twin formula for Fix_SL(sigma) is INCOMPLETE at n=8.
  There must be additional non-twin self-flips at sigma-fixed arcs.

  The first cycle type with TWO 3-cycles is [3,3,1,1] at n=8.
  Non-twin self-flips may arise from automorphisms that exchange
  the two 3-cycles while swapping the fixed points.
""")

# Compute the correction needed at n=8
print("  CORRECTION ANALYSIS for n=8:")
T_8 = 188288
SL_twin_8 = 3145
deficit = (T_8 - SL_twin_8) % 2
print(f"    T_8 = {T_8}, SL_twin = {SL_twin_8}")
print(f"    T - SL_twin = {T_8 - SL_twin_8} (parity: {'odd' if deficit else 'even'})")
print(f"    Need SL correction of at least {deficit} (must be odd)")
print(f"    True SL_8 >= {SL_twin_8 + deficit}")
