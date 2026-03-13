#!/usr/bin/env python3
"""
ocf_directed_fix.py -- kind-pasteur-2026-03-13-S61

Fix the OCF computation: Omega(T) has one vertex per DIRECTED odd cycle,
not per vertex set. Two directed cycles are adjacent iff they share a
tournament vertex.

Key insight from n=7: H_4+H_6 = 2*c7_dir - 16.5 exactly!
The coefficient 2 comes from: each directed Ham cycle uses all n vertices,
conflicts with everything, so adds 2 to I(Omega, 2).

This script:
1. Counts DIRECTED odd cycles properly
2. Verifies H = I(Omega, 2) with directed cycle vertices
3. Checks the "coefficient 2" pattern at n=5 and n=7

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def enumerate_all_directed_cycles(A, n, max_k=None):
    """Enumerate ALL directed odd cycles as ordered tuples (up to rotation).

    Returns list of (vertex_tuple, vertex_set, length) where vertex_tuple
    is the canonical rotation (smallest start) of the directed cycle.

    A directed k-cycle through vertices v_0 -> v_1 -> ... -> v_{k-1} -> v_0
    is represented as the tuple starting with the smallest vertex.
    """
    if max_k is None:
        max_k = n

    cycles = []
    for k in range(3, max_k + 1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            v_k = len(verts)

            # Find ALL directed Hamiltonian cycles on this subset
            # Use DP but track the actual paths
            # dp[mask][v] = list of paths from verts[0] to verts[v] through mask
            # Too expensive for large k — instead just count

            # For counting: use standard DP
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << v_k):
                for v in range(v_k):
                    if not (mask & (1 << v)):
                        continue
                    key = (mask, v)
                    if key not in dp or dp[key] == 0:
                        continue
                    for w in range(v_k):
                        if mask & (1 << w):
                            continue
                        if A[verts[v]][verts[w]]:
                            nk = (mask | (1 << w), w)
                            dp[nk] = dp.get(nk, 0) + dp[key]

            full = (1 << v_k) - 1
            # Cycles returning to verts[0]
            cycle_count_from_0 = 0
            for v in range(1, v_k):
                if (full, v) in dp and A[verts[v]][verts[0]]:
                    cycle_count_from_0 += dp[(full, v)]

            # Each directed cycle on this vertex set that passes through verts[0]
            # is counted once by our DP. But there might be directed cycles on
            # this vertex set that DON'T pass through verts[0]... wait, they all
            # do, since it's a Hamiltonian cycle on the subset.
            # So cycle_count_from_0 = total directed Ham cycles on this subset.
            # Wait, no! Each directed Ham cycle visits ALL vertices including verts[0].
            # So every directed cycle starts from verts[0] in exactly one rotation.
            # But our DP counts directed paths from verts[0] -> ... -> verts[v] -> verts[0].
            # Each directed Hamiltonian cycle on k vertices gives k different starting points,
            # but our DP only uses verts[0] as start. Each directed cycle is counted ONCE.

            # So cycle_count_from_0 = number of directed Hamiltonian cycles on this subset.

            for _ in range(cycle_count_from_0):
                cycles.append((frozenset(subset), k))

    return cycles


def verify_ocf(A, n, max_k=None):
    """Verify H(T) = I(Omega(T), 2) using directed cycles as Omega vertices."""
    H = count_ham_paths(A, n)

    # Get all directed odd cycles (each counted once per vertex set, multiplied)
    cycles = enumerate_all_directed_cycles(A, n, max_k=max_k)
    nc = len(cycles)

    # Build Omega adjacency: two directed cycles are adjacent iff they share a vertex
    # Note: two different directed cycles on the SAME vertex set always share all vertices
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][0] & cycles[j][0]:  # vertex sets intersect
                adj[i][j] = 1
                adj[j][i] = 1

    # Count independent sets
    if nc <= 24:
        nbr = [0] * nc
        for i in range(nc):
            for j in range(nc):
                if adj[i][j]:
                    nbr[i] |= (1 << j)

        alpha = [0] * (nc + 1)
        def backtrack(v, mask, size):
            alpha[size] += 1
            for w in range(v + 1, nc):
                if not (mask & (1 << w)):
                    new_mask = mask | nbr[w]
                    backtrack(w, new_mask, size + 1)

        backtrack(-1, 0, 0)

        I_val = sum(alpha[j] * (2**j) for j in range(nc + 1))
        max_j = max(j for j in range(nc + 1) if alpha[j] > 0)

        return H, I_val, nc, alpha[:max_j+1]
    else:
        # Too large for exact computation
        return H, None, nc, None


# ========================================================================
# ANALYSIS 1: VERIFY OCF AT n=3,4,5 (EXHAUSTIVE)
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: OCF VERIFICATION WITH DIRECTED CYCLES")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    total = 1 << m

    if total > 2**15:
        print(f"\nn={n}: (sampling)")
        import random
        random.seed(42)
        sample = [random.randint(0, total-1) for _ in range(1000)]
    else:
        sample = range(total)

    print(f"\nn={n}: Verifying H = I(Omega_directed, 2)")

    match_count = 0
    mismatch_count = 0
    alpha_by_H = defaultdict(list)

    for bits in sample:
        A = binary_to_tournament(bits, n)
        H, I_val, nc, alpha = verify_ocf(A, n)

        if I_val is not None:
            if H == I_val:
                match_count += 1
            else:
                mismatch_count += 1
                if mismatch_count <= 3:
                    print(f"  MISMATCH: bits={bits}, H={H}, I(Omega,2)={I_val}, nc={nc}")
                    if alpha:
                        print(f"    alpha = {alpha}")

            alpha_by_H[H].append((nc, alpha))

    if mismatch_count == 0:
        print(f"  ALL {match_count} tournaments MATCH! H = I(Omega_directed, 2)")
    else:
        print(f"  {match_count} match, {mismatch_count} MISMATCH")

    # Show alpha decomposition for each H value
    print(f"\n  Alpha decomposition by H:")
    for H in sorted(alpha_by_H.keys()):
        group = alpha_by_H[H]
        nc_set = set(nc for nc, _ in group)
        alpha_example = group[0][1]
        I_check = sum(alpha_example[j] * (2**j) for j in range(len(alpha_example)))
        print(f"    H={H}: nc={sorted(nc_set)}, alpha={alpha_example}, I={I_check}")


# ========================================================================
# ANALYSIS 2: DIRECTED CYCLE COUNT BREAKDOWN AT n=5
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: DIRECTED CYCLE BREAKDOWN AT n=5")
print("=" * 70)

n = 5
m = n * (n - 1) // 2
total = 1 << m

by_score = defaultdict(list)
for bits in range(total):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))

    cycles = enumerate_all_directed_cycles(A, n)
    nc = len(cycles)

    # Count by length
    c3_count = sum(1 for _, k in cycles if k == 3)
    c5_count = sum(1 for _, k in cycles if k == 5)

    by_score[scores].append({
        'bits': bits, 'H': H, 'nc': nc,
        'c3_dir': c3_count, 'c5_dir': c5_count
    })

for sc in sorted(by_score.keys()):
    group = by_score[sc]
    Hs = set(d['H'] for d in group)
    if len(Hs) <= 1 and len(group) > 5:
        # Only show varying classes and small ones
        d = group[0]
        print(f"  Score {sc}: H={sorted(Hs)}, nc={d['nc']}, c3_dir={d['c3_dir']}, c5_dir={d['c5_dir']}")
        continue

    print(f"\n  Score {sc}:")
    by_H_local = defaultdict(list)
    for d in group:
        by_H_local[d['H']].append(d)

    for H in sorted(by_H_local.keys()):
        g = by_H_local[H]
        nc_set = set(d['nc'] for d in g)
        c3_set = set(d['c3_dir'] for d in g)
        c5_set = set(d['c5_dir'] for d in g)
        print(f"    H={H}: count={len(g)}, nc={sorted(nc_set)}, "
              f"c3_dir={sorted(c3_set)}, c5_dir={sorted(c5_set)}")


# ========================================================================
# ANALYSIS 3: THE COEFFICIENT-2 THEOREM
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: THE COEFFICIENT-2 THEOREM")
print("=" * 70)

print("""
THEOREM (Coefficient-2 for longest cycles):
  For a directed k-cycle C using ALL n vertices (Hamiltonian cycle, k=n),
  C conflicts with every other cycle in Omega(T) since it uses all vertices.
  Therefore, adding C to Omega increases alpha_1 by 1 but leaves alpha_j
  unchanged for j >= 2.

  Consequence: I(Omega, 2) increases by exactly 2^1 = 2.

  This explains:
  - At n=5, score class (1,2,2,2,3): H = 9 + 2*c5_dir
    (c5 are Hamiltonian 5-cycles, each adding 2 to H)
  - At n=7, regular: H = 141 + 2*c7_dir
    (c7 are Hamiltonian 7-cycles, each adding 2 to H)

  The pattern: within a fixed score class (or regular class), the ONLY
  cycle length that varies independently is the HAMILTONIAN one (length n).
  All shorter cycle counts are either score-determined (c3) or constrained
  by the graph structure to not affect H independently (c5 at n=7).
""")

# Verify at n=5 explicitly
n = 5
m = n * (n - 1) // 2
total = 1 << m

print(f"Verification at n={n}, score class (1,2,2,2,3):")
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = tuple(sorted([sum(A[v]) for v in range(n)]))
    if scores != (1, 2, 2, 2, 3):
        continue

    H = count_ham_paths(A, n)
    cycles = enumerate_all_directed_cycles(A, n)
    nc = len(cycles)
    c5_dir = sum(1 for _, k in cycles if k == 5)

    predicted_H_part = 2 * c5_dir
    residual = H - predicted_H_part

    # Print first few of each H
    if bits < 200 or H in [11, 13, 15]:
        print(f"  bits={bits}: H={H}, c5_dir={c5_dir}, 2*c5_dir={predicted_H_part}, "
              f"H - 2*c5_dir = {residual}")


# ========================================================================
# ANALYSIS 4: REGULAR n=7 — VERIFY WITH DIRECTED CYCLES
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: REGULAR n=7 WITH DIRECTED CYCLE COUNTING")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
total = 1 << m

print(f"Scanning {total} tournaments for regular ones...")

by_H = defaultdict(list)
checked = 0
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if any(s != 3 for s in scores):
        continue

    H = count_ham_paths(A, n)

    # Count directed cycles
    cycles = enumerate_all_directed_cycles(A, n)
    nc = len(cycles)
    c3_dir = sum(1 for _, k in cycles if k == 3)
    c5_dir = sum(1 for _, k in cycles if k == 5)
    c7_dir = sum(1 for _, k in cycles if k == 7)

    by_H[H].append({
        'nc': nc, 'c3_dir': c3_dir, 'c5_dir': c5_dir, 'c7_dir': c7_dir
    })
    checked += 1

print(f"Found {checked} regular tournaments")
print(f"\nDirected cycle decomposition:")
print(f"  {'H':>5s} | {'cnt':>5s} | {'nc':>5s} | {'c3_dir':>7s} | {'c5_dir':>7s} | {'c7_dir':>7s} | {'H-2*c7':>8s}")
print(f"  {'':->5s}-+-{'':->5s}-+-{'':->5s}-+-{'':->7s}-+-{'':->7s}-+-{'':->7s}-+-{'':->8s}")

for H in sorted(by_H.keys()):
    g = by_H[H]
    def show(key):
        vals = sorted(set(d[key] for d in g))
        return str(vals[0]) if len(vals) == 1 else str(vals)

    c7_vals = list(set(d['c7_dir'] for d in g))
    residual = H - 2*c7_vals[0] if len(c7_vals) == 1 else "?"

    print(f"  {H:>5d} | {len(g):>5d} | {show('nc'):>5s} | {show('c3_dir'):>7s} | "
          f"{show('c5_dir'):>7s} | {show('c7_dir'):>7s} | {residual}")

# Verify OCF for one representative of each class
print(f"\nOCF verification (one per class):")
for H in sorted(by_H.keys()):
    d = by_H[H][0]
    nc = d['nc']

    # We need the actual tournament for OCF verification
    # Find the first tournament with this H
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        scores = [sum(A[v]) for v in range(n)]
        if any(s != 3 for s in scores):
            continue
        if count_ham_paths(A, n) != H:
            continue

        # Full OCF verification
        H_actual, I_val, nc_actual, alpha = verify_ocf(A, n)
        if I_val is not None:
            print(f"  H={H}: I(Omega_dir, 2) = {I_val}, nc={nc_actual}, "
                  f"alpha={alpha}, match={'YES' if H==I_val else 'NO'}")
        else:
            print(f"  H={H}: nc={nc_actual} too large for exact IS enumeration")
        break


print("\n" + "=" * 70)
print("DONE.")
