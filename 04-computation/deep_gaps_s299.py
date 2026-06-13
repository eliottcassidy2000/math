#!/usr/bin/env python3
"""
deep_gaps_s299.py — The 129 order-≥3 gaps: what fills them?
opus-2026-03-24-S299

At n=6: 129 class pairs need order ≥ 3 tile flips.
Compute waggly-3 (order 3) to see how many it fills.
Then check if 1+2+3 = complete, or if order 4+ is needed.
Look for recursive/structural patterns in the gaps.
"""

import sys, time
from math import comb
from itertools import combinations
from collections import defaultdict, Counter
import numpy as np

# Reuse build_data from earlier
from wiggly_gaps_s296 import build_data

sys.stdout.reconfigure(line_buffering=True)

def analyze_deep_gaps(n):
    print(f"\n{'#'*72}")
    print(f"  n = {n}: DEEP GAPS (ORDER ≥ 3)")
    print(f"{'#'*72}")

    t0 = time.time()
    V, mi, tc, m, total = build_data(n)
    print(f"  V={V}, m={m}, 2^m={total}")

    # Build edge sets for orders 1, 2, 3
    E = [set() for _ in range(m+1)]  # E[k] = edges from order-k flips

    for tb in range(total):
        ma = tc[tb]
        if ma < 0: continue

        # Order 1 (wiggly)
        for ti in range(m):
            mb = tc[tb ^ (1 << (m-1-ti))]
            if mb >= 0 and ma != mb:
                E[1].add((min(ma,mb), max(ma,mb)))

        # Order 2 (waggly-2)
        for ti in range(m):
            for tj in range(ti+1, m):
                mb = tc[tb ^ (1 << (m-1-ti)) ^ (1 << (m-1-tj))]
                if mb >= 0 and ma != mb:
                    E[2].add((min(ma,mb), max(ma,mb)))

        # Order 3 (waggly-3)
        for ti in range(m):
            for tj in range(ti+1, m):
                for tk in range(tj+1, m):
                    mb = tc[tb ^ (1 << (m-1-ti)) ^ (1 << (m-1-tj)) ^ (1 << (m-1-tk))]
                    if mb >= 0 and ma != mb:
                        E[3].add((min(ma,mb), max(ma,mb)))

        # Order m (complement)
        mb = tc[((1 << m) - 1) ^ tb]
        if mb >= 0 and ma != mb:
            E[m].add((min(ma,mb), max(ma,mb)))

    total_pairs = V*(V-1)//2
    dt = time.time() - t0

    print(f"  Computed in {dt:.1f}s")

    # Cumulative filling
    F = [set()]  # F[0] = empty
    cumulative = set()
    for k in range(1, min(m+1, 6)):
        if k <= 3 or k == m:
            cumulative = cumulative | E[k]
            F.append(cumulative.copy())
            remaining = total_pairs - len(cumulative)
            print(f"  F_{k} (orders 1..{k}): {len(cumulative)}/{total_pairs} edges "
                  f"({100*len(cumulative)/total_pairs:.1f}%), {remaining} gaps remain")

    # Also add complement
    all_with_comp = cumulative | E[m]
    remaining_with_comp = total_pairs - len(all_with_comp)
    print(f"  F_3 + complement: {len(all_with_comp)}/{total_pairs} "
          f"({100*len(all_with_comp)/total_pairs:.1f}%), {remaining_with_comp} gaps remain")

    # The ORDER-3-ONLY edges (not in order 1 or 2)
    E3_only = E[3] - E[1] - E[2]
    E3_new = E[3] - E[1] - E[2]  # same thing

    print(f"\n  ORDER-3 ANALYSIS:")
    print(f"    E[3] total: {len(E[3])} edges")
    print(f"    E[3] new (not in 1 or 2): {len(E3_new)}")
    print(f"    These fill the distance-3 gaps from S296")

    # Gaps after order 3
    gaps_after_3 = total_pairs - len(E[1] | E[2] | E[3])
    print(f"\n    Gaps after orders 1+2+3: {gaps_after_3}")

    if gaps_after_3 == 0:
        print(f"    *** ORDERS 1+2+3 ARE COMPLETE! ***")
        print(f"    Every class pair is reachable by ≤ 3 tile flips.")
    else:
        print(f"    {gaps_after_3} pairs still need order ≥ 4.")

    # Analyze the remaining gaps (if any)
    if gaps_after_3 > 0:
        remaining_gaps = []
        all_covered = E[1] | E[2] | E[3]
        for i in range(V):
            for j in range(i+1, V):
                if (i,j) not in all_covered:
                    remaining_gaps.append((i, j))

        print(f"\n    ORDER-≥4 GAPS:")
        for (i, j) in remaining_gaps[:20]:
            print(f"      ({i},{j}): H=({mi[i]['H']},{mi[j]['H']}), "
                  f"SC=({'SC' if mi[i]['sc'] else 'NS'},{'SC' if mi[j]['sc'] else 'NS'})")

        # Check: are these filled by order 4?
        if n <= 6 and gaps_after_3 <= 50:
            print(f"\n    Checking order 4...")
            E4_fills = set()
            for tb in range(total):
                ma = tc[tb]
                if ma < 0: continue
                for combo in combinations(range(m), 4):
                    flip_mask = 0
                    for ti in combo:
                        flip_mask ^= (1 << (m-1-ti))
                    mb = tc[tb ^ flip_mask]
                    if mb >= 0 and ma != mb:
                        edge = (min(ma,mb), max(ma,mb))
                        if edge not in all_covered:
                            E4_fills.add(edge)

            print(f"    Order-4 fills: {len(E4_fills)} of {gaps_after_3} remaining gaps")
            gaps_after_4 = gaps_after_3 - len(E4_fills)
            print(f"    Gaps after 1+2+3+4: {gaps_after_4}")

            if gaps_after_4 == 0:
                print(f"    *** ORDERS 1+2+3+4 ARE COMPLETE! ***")

    # ΔH analysis by order
    print(f"\n  ΔH BY ORDER:")
    for k in [1, 2, 3]:
        if E[k]:
            new_at_k = E[k] - (E[1] | E[2] | E[3] if k > 1 else set())
            # Actually compute new edges AT each order
            if k == 1:
                new_edges = E[1]
            elif k == 2:
                new_edges = E[2] - E[1]
            else:
                new_edges = E[3] - E[1] - E[2]

            if new_edges:
                dH = [abs(mi[a]['H'] - mi[b]['H']) for (a,b) in new_edges]
                print(f"    Order {k} new edges: {len(new_edges)}, "
                      f"mean |ΔH| = {np.mean(dH):.1f}, range [{min(dH)}, {max(dH)}]")

    # RECURSIVE STRUCTURE: do gaps correspond to Mode B changes?
    # Mode B: inner sub-tournament on vertices {1,...,n-2}
    # If two classes differ in their inner structure by ≥ 3 tiles,
    # they're in the order-≥3 gap.
    print(f"\n  STRUCTURAL ANALYSIS OF GAPS:")

    # Compute min Hamming distances between all gap pairs
    class_tilings = defaultdict(list)
    for tb in range(total):
        if tc[tb] >= 0:
            class_tilings[tc[tb]].append(tb)

    gap_pairs = [(i,j) for i in range(V) for j in range(i+1,V)
                 if (i,j) not in (E[1] | E[2])]

    if len(gap_pairs) <= 200:
        dist_3_gaps = []
        dist_4plus_gaps = []

        for (i, j) in gap_pairs:
            min_hd = m + 1
            for ti in class_tilings[i]:
                for tj in class_tilings[j]:
                    hd = bin(ti ^ tj).count('1')
                    min_hd = min(min_hd, hd)
                    if min_hd <= 2: break
                if min_hd <= 2: break

            if min_hd == 3:
                dist_3_gaps.append((i, j))
            elif min_hd >= 4:
                dist_4plus_gaps.append((i, j, min_hd))

        print(f"    Distance-3 gaps: {len(dist_3_gaps)}")
        print(f"    Distance-4+ gaps: {len(dist_4plus_gaps)}")
        if dist_4plus_gaps:
            print(f"    Distance-4+ details:")
            for (i, j, d) in dist_4plus_gaps:
                print(f"      ({i},{j}): distance={d}, H=({mi[i]['H']},{mi[j]['H']})")

    return V

if __name__ == '__main__':
    print("=" * 72)
    print("  DEEP GAPS: WHAT FILLS ORDER ≥ 3?")
    print("  opus-2026-03-24-S299")
    print("=" * 72)

    # Need to be in the right directory
    import os
    os.chdir('/Users/e/Documents/GitHub/math')
    sys.path.insert(0, '04-computation')

    for n in [5, 6]:
        analyze_deep_gaps(n)

    print(f"\n{'='*72}")
    print("  THE COMPLETENESS THRESHOLD")
    print(f"{'='*72}")
    print("""
  The COMPLETENESS ORDER k* = minimum k such that
  F_k = orders 1..k fills ALL class pairs.

  n=4: k* = 1 (wiggly alone is complete)
  n=5: k* = ?
  n=6: k* = ?

  If k* ≤ diameter of the metagraph, then every gap
  is accountable by the Hamming distance structure.

  RECURSIVE STRUCTURE:
  The gaps may correspond to classes differing in their
  MODE B inner structure by a coordinated multi-tile change.
  This would connect the gap order to the recursive depth
  of the staircase decomposition.
""")

    print("DONE.")
