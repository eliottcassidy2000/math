#!/usr/bin/env python3
"""
sl_mine_bruteforce_s20eh.py — Brute-force R(ct) for EVERY cycle type at n=3..7
kind-pasteur-2026-03-23-S20eh

For each cycle type ct with at least one even cycle:
  R(ct) = #{labeled (T, e) : sigma(T flip e) = T}
where sigma is a representative permutation with cycle type ct.

This gives SL_mine(n) = (1/n!) * sum_ct count(ct) * R(ct).

The key insight from the EXACT analysis:
- 2-cycle (u,v): arc {u,v} is a FIXED POINT of sigma on arcs, but sigma
  REVERSES its direction. So it ALWAYS has defect 1.
- k-cycle for k>=3: arcs within the cycle form orbits of various sizes.
  The defect depends on the orbit structure.
- For defect EXACTLY 1: we need sigma(T flip e) = T, i.e., flipping arc e
  in T and then applying sigma gives T back.

This is equivalent to: sigma(T) differs from T at exactly one arc position,
AND that position is sigma(e), AND flipping e compensates.
"""

import sys
from math import factorial, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  BRUTE-FORCE R(ct) PER CYCLE TYPE")
print("  kind-pasteur-2026-03-23-S20eh")
print("=" * 80)

def partitions(n, max_val=None):
    if max_val is None: max_val = n
    if n == 0: yield (); return
    for i in range(min(n, max_val), 0, -1):
        for rest in partitions(n - i, i):
            yield (i,) + rest

def count_perms(ct, n):
    mult = Counter(ct)
    d = 1
    for c, m in mult.items(): d *= (c**m) * factorial(m)
    return factorial(n) // d

def build_sigma(ct, n):
    """Build representative permutation with cycle type ct."""
    sigma = [0] * n
    pos = 0
    for c in ct:
        for i in range(c - 1):
            sigma[pos + i] = pos + i + 1
        sigma[pos + c - 1] = pos
        pos += c
    return sigma

def apply_sigma_to_tournament(bits, sigma, n, ALL_ARCS, arc_idx):
    """Apply sigma to tournament T (as bit vector). Returns sigma(T) bits."""
    new_bits = 0
    for k, (i, j) in enumerate(ALL_ARCS):
        # Arc (i,j): T has i->j if bit k set, else j->i
        si, sj = sigma[i], sigma[j]
        # sigma maps arc direction: if T has i->j, then sigma(T) has si->sj
        if bits & (1 << k):
            # i -> j in T, so si -> sj in sigma(T)
            if si < sj:
                target_k = arc_idx[(si, sj)]
                new_bits |= (1 << target_k)
            else:
                # si > sj: the arc position is (sj, si), direction is si->sj
                # which means sj is NOT -> si, so bit for (sj,si) is 0
                pass  # bit stays 0
        else:
            # j -> i in T, so sj -> si in sigma(T)
            if sj < si:
                target_k = arc_idx[(sj, si)]
                new_bits |= (1 << target_k)
            else:
                pass
    return new_bits

def compute_R_bruteforce(ct, n):
    """Brute-force: count (T, e) pairs where sigma(T flip e) = T."""
    sigma = build_sigma(ct, n)
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    arc_idx = {a: i for i, a in enumerate(ALL_ARCS)}

    count = 0
    for bits in range(1 << m):
        for k in range(m):
            flipped = bits ^ (1 << k)
            sigma_flipped = apply_sigma_to_tournament(flipped, sigma, n, ALL_ARCS, arc_idx)
            if sigma_flipped == bits:
                count += 1
    return count


# Known values
T_known = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
E_known = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

for n in range(3, 8):
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)

    print(f"\n{'='*60}")
    print(f"  n = {n}, m = {m} arcs, 2^m = {2**m} tournaments")
    print(f"{'='*60}")

    total_raw = 0
    ct_details = []

    for ct in partitions(n):
        has_even = any(c % 2 == 0 for c in ct)
        if not has_even:
            continue

        nperms = count_perms(ct, n)
        R = compute_R_bruteforce(ct, n)
        contrib = nperms * R
        total_raw += contrib

        ct_details.append((ct, R, nperms, contrib))
        if R > 0:
            print(f"  ct={ct}: R={R}, perms={nperms}, total={contrib}")

    sl_mine = total_raw // factorial(n)
    sl_expected = T_known.get(n, 0) - 2 * E_known.get(n, 0) if n in T_known and n in E_known else None

    print(f"\n  Total raw = {total_raw}")
    print(f"  SL_mine = {total_raw} / {factorial(n)} = {sl_mine}")
    if sl_expected is not None:
        print(f"  Expected SL_mine = T - 2E = {T_known[n]} - 2*{E_known[n]} = {sl_expected}")
        print(f"  {'MATCH!' if sl_mine == sl_expected else f'OFF BY {sl_mine - sl_expected}'}")

    # Show ALL cycle types (including zero contribution)
    print(f"\n  ALL cycle types with even parts:")
    for ct, R, nperms, contrib in sorted(ct_details, key=lambda x: -x[2]):
        even_parts = [c for c in ct if c % 2 == 0]
        print(f"    {str(ct):20s}  even={even_parts}  R={R:6d}  perms={nperms:6d}  contrib={contrib:8d}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")

    if n >= 7 and elapsed > 60:
        print("  (Skipping larger n due to time)")
        break

print("\nDONE.")
print("=" * 80)
