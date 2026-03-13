#!/usr/bin/env python3
"""
overlap_weight_counts.py -- kind-pasteur-2026-03-13-S60

FAST VERSION: Just compute overlap weight distribution by cycle count per orbit.

At p=11 Paley, we observed 5-cycle orbits have:
  count=1: w=21152 or 21158
  count=2: w=21155 or 21161
  count=3: w=21156

The OVERLAP WEIGHT of a 5-cycle orbit depends on the number of directed
cycles it supports. Does w = a + b*n(V) for some constants a, b?

More precisely: for a 5-element subset V supporting n(V) directed cycles,
the total overlap weight of each directed cycle through V is w.
Is w a function of n(V) alone?
"""

from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def fast_overlap_analysis(p, S, name):
    """Compute overlap weights grouped by n(V) per vertex set."""
    A = build_adj(p, S)

    # For each vertex set of size k, compute n(V) and overlap data
    for k in range(3, p + 1, 2):
        # Build list of (vertex_set, n(V))
        vset_data = []
        for subset in combinations(range(p), k):
            fs = frozenset(subset)
            nc = count_ham_cycles(A, list(subset))
            if nc > 0:
                vset_data.append((fs, nc))

        if not vset_data:
            continue

        c_k_total = sum(nc for _, nc in vset_data)

        # For each active vertex set, compute its "overlap" with all other
        # active vertex sets
        active_sets = [(fs, nc) for fs, nc in vset_data]

        # For each active set V, compute:
        # - How many other active sets share at least one vertex
        # - Total n(V') for overlapping sets (total conflicting cycles)
        overlaps = []
        for i, (fs1, nc1) in enumerate(active_sets):
            n_conflicting_vsets = 0
            total_conflicting_cycles = 0
            for j, (fs2, nc2) in enumerate(active_sets):
                if i == j:
                    continue
                if fs1 & fs2:
                    n_conflicting_vsets += 1
                    total_conflicting_cycles += nc2

            overlaps.append({
                'vset': fs1,
                'n': nc1,
                'conflicting_vsets': n_conflicting_vsets,
                'conflicting_cycles': total_conflicting_cycles,
            })

        # Group by n(V)
        by_n = defaultdict(list)
        for ov in overlaps:
            by_n[ov['n']].append(ov)

        print(f"\n  k={k}: {len(active_sets)} active sets, c_{k}={c_k_total}")
        print(f"  {'n(V)':>5} {'count':>6} {'conf_vsets':>12} {'conf_cycles':>13} {'w_per_cycle':>13}")

        for n in sorted(by_n.keys()):
            group = by_n[n]
            cvs = set(ov['conflicting_vsets'] for ov in group)
            ccs = set(ov['conflicting_cycles'] for ov in group)

            # w_per_cycle: total conflicting directed cycles
            # This is the overlap weight per directed cycle FROM this vertex set
            # If n(V)=n, each directed cycle on V conflicts with total_conflicting_cycles others
            # PLUS (n-1) other directed cycles on V itself
            w_full = set(ov['conflicting_cycles'] + (ov['n'] - 1) for ov in group)

            if len(cvs) == 1 and len(ccs) == 1:
                ov = group[0]
                print(f"  {n:>5} {len(group):>6} {ov['conflicting_vsets']:>12} "
                      f"{ov['conflicting_cycles']:>13} {list(w_full)[0]:>13}  CONSTANT")
            else:
                min_cc = min(ov['conflicting_cycles'] for ov in group)
                max_cc = max(ov['conflicting_cycles'] for ov in group)
                print(f"  {n:>5} {len(group):>6} {'['+str(min(cvs))+'-'+str(max(cvs))+']':>12} "
                      f"{'['+str(min_cc)+'-'+str(max_cc)+']':>13} "
                      f"{'['+str(min(w_full))+'-'+str(max(w_full))+']':>13}  VARIES")


for p_val in [7, 11]:
    m = (p_val - 1) // 2
    S_qr = sorted(j for j in range(1, p_val) if pow(j, (p_val - 1) // 2, p_val) == 1)
    S_int = list(range(1, m + 1))

    print(f"\n{'='*70}")
    print(f"p={p_val}")
    print(f"{'='*70}")

    for name, S in [("Paley", S_qr), ("Interval", S_int)]:
        print(f"\n--- {name} S={S} ---")
        fast_overlap_analysis(p_val, S, name)

print("\nDONE.")
