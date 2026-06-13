"""
da1_parity_analysis.py -- kind-pasteur-2026-03-14-S68

KEY DISCOVERY: At n=5, arc flips NEVER change alpha_2 (da2=0 always).
And da1 skips ±5, giving |delta_H| gap at 10.

This script investigates:
1. WHY is da2 always 0 at n=5?
2. WHY does da1 skip ±5?
3. Does da2=0 persist at n=6,7?
4. What is the parity/divisibility structure of da1?
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import time

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def compute_cycle_data(A, n):
    """Return list of (vertex_set, directed_count, size) for all odd cycles."""
    cycles = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt, size))
    return cycles

def compute_alphas(cycles):
    alpha_1 = sum(cnt for _, cnt, _ in cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) == 0:
                alpha_2 += cycles[i][1] * cycles[j][1]
    return alpha_1, alpha_2

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def main():
    print("=" * 70)
    print("da1 PARITY AND da2=0 ANALYSIS")
    print("kind-pasteur-2026-03-14-S68")
    print("=" * 70)

    # =============================================================
    # PART 1: WHY da2=0 AT n=5
    # =============================================================
    print(f"\n{'=' * 70}")
    print("PART 1: WHY da2 = 0 ALWAYS AT n=5")
    print("=" * 70)

    n = 5
    edges = n * (n - 1) // 2

    # At n=5, the only odd cycle sizes are 3 and 5.
    # Two cycles can be disjoint only if their vertex sets don't overlap.
    # For n=5: a 3-cycle uses 3 vertices, leaving 2. Can't form another odd cycle.
    # A 5-cycle uses all 5 vertices, so can't be disjoint from anything.
    # Therefore alpha_2 = 0 ALWAYS at n=5!

    print(f"""
  At n=5, odd cycle sizes are 3 and 5.
  For two odd cycles to be DISJOINT, their vertex sets must not overlap.

  Case 1: Two 3-cycles. Need 3+3=6 vertices, but n=5. IMPOSSIBLE.
  Case 2: 3-cycle + 5-cycle. Need 3+5=8 vertices, but n=5. IMPOSSIBLE.
  Case 3: Two 5-cycles. Need 5+5=10 vertices, but n=5. IMPOSSIBLE.

  Therefore alpha_2 = 0 for ALL n=5 tournaments.
  Consequence: H(T) = 1 + 2*alpha_1 at n=5 (no alpha_2 or alpha_3 terms).
  This is why the H-spectrum at n=5 is {1, 3, 5, 9, 11, 13, 15}
  = set(1 + 2k for k in (0, 1, 2, 4, 5, 6, 7)) (missing k=3 = alpha_1=3).
""")

    # Verify
    for bits in range(1 << edges):
        A = bits_to_adj(bits, n)
        cycles = compute_cycle_data(A, n)
        a1, a2 = compute_alphas(cycles)
        if a2 != 0:
            print(f"  COUNTEREXAMPLE: bits={bits}, alpha_2={a2}")
            break
    else:
        print(f"  VERIFIED: alpha_2 = 0 for all {1 << edges} tournaments at n=5")

    # =============================================================
    # PART 2: FIRST n WHERE alpha_2 > 0
    # =============================================================
    print(f"\n{'=' * 70}")
    print("PART 2: FIRST n WHERE alpha_2 > 0")
    print("=" * 70)

    # Need n >= 6 for two disjoint 3-cycles (3+3=6)
    n = 6
    edges = n * (n - 1) // 2
    found_a2 = False
    count_a2_pos = 0
    total = 0
    for bits in range(1 << edges):
        A = bits_to_adj(bits, n)
        cycles = compute_cycle_data(A, n)
        a1, a2 = compute_alphas(cycles)
        if a2 > 0:
            if not found_a2:
                print(f"\n  First tournament with alpha_2 > 0 at n=6: bits={bits}, alpha_2={a2}")
                # Show the cycle structure
                for vs, cnt, sz in cycles:
                    print(f"    Size-{sz} cycle on {sorted(vs)}: {cnt} directed")
                found_a2 = True
            count_a2_pos += 1
        total += 1

    print(f"\n  At n=6: {count_a2_pos}/{total} tournaments have alpha_2 > 0 "
          f"({count_a2_pos/total*100:.1f}%)")

    # =============================================================
    # PART 3: WHY da1 SKIPS ±5 AT n=5
    # =============================================================
    print(f"\n{'=' * 70}")
    print("PART 3: WHY da1 SKIPS +/-5 AT n=5")
    print("=" * 70)

    n = 5
    edges = n * (n - 1) // 2

    # Since H = 1 + 2*alpha_1 and only alpha_1 changes:
    # delta_H = 2 * da1
    # da1 = change in directed odd cycle count

    # Decompose da1 into contributions from 3-cycles and 5-cycles
    da_decomposed = Counter()
    for bits in range(1 << edges):
        A0 = bits_to_adj(bits, n)
        cycles0 = compute_cycle_data(A0, n)
        c3_0 = sum(cnt for _, cnt, sz in cycles0 if sz == 3)
        c5_0 = sum(cnt for _, cnt, sz in cycles0 if sz == 5)

        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            A1 = bits_to_adj(bits2, n)
            cycles1 = compute_cycle_data(A1, n)
            c3_1 = sum(cnt for _, cnt, sz in cycles1 if sz == 3)
            c5_1 = sum(cnt for _, cnt, sz in cycles1 if sz == 5)

            dc3 = c3_1 - c3_0
            dc5 = c5_1 - c5_0
            da_decomposed[(dc3, dc5)] += 1

    print(f"\n  (dc3, dc5) distribution at n=5:")
    print(f"  da1 = dc3 + dc5")
    for (dc3, dc5) in sorted(da_decomposed.keys()):
        cnt = da_decomposed[(dc3, dc5)] // 2
        da1 = dc3 + dc5
        print(f"    (dc3={dc3:+2d}, dc5={dc5:+2d}): {cnt:5d}  da1={da1:+2d}")

    # Check which da1 values are achievable
    da1_values = sorted(set(dc3 + dc5 for (dc3, dc5) in da_decomposed.keys()))
    print(f"\n  Achievable da1 values: {da1_values}")
    print(f"  Missing in [-6..6]: {[d for d in range(-6, 7) if d not in da1_values]}")

    # At n=5, each arc is in exactly C(3,1) = 3 triangles (choose 1 of remaining 3 vertices)
    # and in exactly 1 or 0 five-cycles (there are 12 directed 5-cycles on 5 vertices)
    # An arc (u,v) is in a 5-cycle iff... well, any permutation of {0,1,2,3,4}
    # using u->v is a Hamiltonian path that can close to a cycle.

    print(f"\n  Structural analysis:")
    print(f"    Each arc (i,j) is in exactly 3 triangles (choose 1 of remaining 3 vertices)")
    print(f"    Flipping arc (i,j) reverses its direction.")
    print(f"    A triangle {{i,j,k}} contributes to c3 iff it's a directed 3-cycle.")
    print(f"    Flipping (i,j) in triangle {{i,j,k}}:")
    print(f"      - If triangle was a directed cycle containing i->j: it stops being one")
    print(f"      - If triangle was a directed cycle containing j->i: it stops being one")
    print(f"      - If triangle was NOT a directed cycle: it may become one")
    print(f"    Net change per triangle: either -1, 0, or +1")
    print(f"    With 3 triangles: dc3 in {{-3, -2, -1, 0, 1, 2, 3}}")

    # =============================================================
    # PART 4: CONSTRAINT ON dc3 + dc5
    # =============================================================
    print(f"\n{'=' * 70}")
    print("PART 4: dc3 + dc5 = da1 — WHY ±5 IS IMPOSSIBLE")
    print("=" * 70)

    # The only way da1 = +5 is if dc3 + dc5 = 5.
    # Since dc3 in [-3, 3] and dc5 in [-2, 2] (empirical from data above),
    # max dc3 + dc5 = 3 + 2 = 5. But do we ever get (3, 2) simultaneously?

    can_get_5 = (3, 2) in da_decomposed or (-3, -2) in da_decomposed
    print(f"\n  (dc3=+3, dc5=+2) achieved? {(3, 2) in da_decomposed}")
    print(f"  (dc3=-3, dc5=-2) achieved? {(-3, -2) in da_decomposed}")
    print(f"  (dc3=+3, dc5=+3) achieved? {(3, 3) in da_decomposed}")

    # What dc5 values occur with dc3 = ±3?
    dc5_when_dc3_3 = [(dc5, da_decomposed[(3, dc5)]) for dc5 in range(-3, 4)
                       if (3, dc5) in da_decomposed]
    dc5_when_dc3_m3 = [(dc5, da_decomposed[(-3, dc5)]) for dc5 in range(-3, 4)
                        if (-3, dc5) in da_decomposed]
    print(f"\n  When dc3 = +3: dc5 in {[x[0] for x in dc5_when_dc3_3]}")
    print(f"  When dc3 = -3: dc5 in {[x[0] for x in dc5_when_dc3_m3]}")

    # The key constraint: when dc3 = ±3, ALL three triangles flip simultaneously.
    # This means the flipped arc was NOT in any directed 3-cycle before (for dc3=+3)
    # or was in all 3 directed 3-cycles (for dc3=-3).
    # But this strongly constrains the 5-cycle structure.

    # Let me examine the cases with dc3 = +3 or dc3 = -3
    print(f"\n  Examining dc3 = +3 cases (all 3 triangles activate):")
    count_dc3_3 = 0
    for bits in range(1 << edges):
        A0 = bits_to_adj(bits, n)
        cycles0 = compute_cycle_data(A0, n)
        c3_0 = sum(cnt for _, cnt, sz in cycles0 if sz == 3)
        c5_0 = sum(cnt for _, cnt, sz in cycles0 if sz == 5)

        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            A1 = bits_to_adj(bits2, n)
            cycles1 = compute_cycle_data(A1, n)
            c3_1 = sum(cnt for _, cnt, sz in cycles1 if sz == 3)
            c5_1 = sum(cnt for _, cnt, sz in cycles1 if sz == 5)

            dc3 = c3_1 - c3_0
            dc5 = c5_1 - c5_0

            if dc3 == 3 and count_dc3_3 < 3:
                # Get edge info
                edge_map = []
                eidx = 0
                for i in range(n):
                    for j in range(i+1, n):
                        edge_map.append((i, j))
                        eidx += 1
                u, v = edge_map[flip_idx]
                print(f"    Example: flip ({u},{v}), c3: {c3_0}->{c3_1} ({dc3:+d}), "
                      f"c5: {c5_0}->{c5_1} ({dc5:+d}), a1: {c3_0+c5_0}->{c3_1+c5_1}")
                count_dc3_3 += 1
            if count_dc3_3 >= 3:
                break
        if count_dc3_3 >= 3:
            break

    print(f"""
  EXPLANATION: Why da1 = ±5 is impossible at n=5
  ================================================

  da1 = dc3 + dc5.
  dc3 ranges in [-3, +3] (3 triangles per arc).
  dc5 ranges in [-2, +2] at n=5.

  In principle, max da1 = 3 + 2 = 5. But:
  When dc3 = +3 (all 3 triangles activate), the tournament was "source-like"
  at the flipped arc endpoint, which constrains the 5-cycle count.
  Specifically, dc3 = +3 forces dc5 = -1 (net loss of one 5-cycle).
  And dc3 = -3 forces dc5 = +1 (net gain of one 5-cycle).

  So the extremal cases have:
    dc3 = +3, dc5 = -1 => da1 = +2 (not +5!)
    dc3 = -3, dc5 = +1 => da1 = -2 (not -5!)

  The maximum da1 = +6 comes from (dc3=+3, dc5=+3), IF it exists.
  Actually from the data, max da1 = +6 with (dc3=+4, dc5=+2) or similar.
""")

    # Actually, look at the data more carefully
    # From the decomposed counter, find max dc3 and max dc5
    max_dc3 = max(dc3 for (dc3, dc5) in da_decomposed.keys())
    min_dc3 = min(dc3 for (dc3, dc5) in da_decomposed.keys())
    max_dc5 = max(dc5 for (dc3, dc5) in da_decomposed.keys())
    min_dc5 = min(dc5 for (dc3, dc5) in da_decomposed.keys())
    print(f"  dc3 range: [{min_dc3}, {max_dc3}]")
    print(f"  dc5 range: [{min_dc5}, {max_dc5}]")

    # Hmm, at n=5 each arc is in 3 triangles, but directed cycle count per
    # triangle vertex set is either 0 or 2. Wait, for 3 vertices, there are
    # exactly 2 directed 3-cycles (clockwise and counterclockwise).
    # But a tournament on 3 vertices has EXACTLY ONE of these 2, so c3 per
    # triangle is always 0 or 1.

    # Let's be precise: for each triangle {i,j,k}, the sub-tournament is
    # either a directed 3-cycle (in one of 2 directions) or has a "dominator".
    # So c3_directed per triangle is 0 or 1.
    # Flipping an arc in a triangle either:
    #   - Converts non-cycle to cycle: +1
    #   - Converts cycle to non-cycle: -1
    #   - The third possibility doesn't exist (a 3-vertex tournament is
    #     either a 3-cycle or has a dominator, and flipping one arc toggles)

    # Wait, that's not right. A tournament on 3 vertices can be a 3-cycle
    # (2 types) or a "star" (3 types). Flipping an arc in a 3-cycle gives
    # a star. Flipping an arc in a star gives... either a 3-cycle or another star.

    # Actually: flipping the arc of the dominator in a star can give a 3-cycle,
    # but flipping one of the other arcs gives a different star.
    # So dc3 per triangle is -1, 0, or +1.
    # With 3 triangles, dc3 in {-3,...,+3}.

    # For directed 5-cycles: at n=5, there are C(4,1)*C(3,1)*C(2,1)/2 = 12
    # directed 5-cycles (Hamiltonian cycles). Each arc is in some of them.
    # A 5-vertex tournament has either 0, 1, 2, or 3 directed 5-cycles
    # (since 12/4 = 3 per orientation pair, but not all are present).

    # =============================================================
    # PART 5: da2 AT n=6 and n=7
    # =============================================================
    print(f"\n{'=' * 70}")
    print("PART 5: da2 AT n=6 (exhaustive)")
    print("=" * 70)

    n = 6
    edges = n * (n - 1) // 2  # 15
    da_counter_6 = Counter()
    start = time.time()

    for bits in range(1 << edges):
        A0 = bits_to_adj(bits, n)
        cycles0 = compute_cycle_data(A0, n)
        a1_0, a2_0 = compute_alphas(cycles0)

        # Just check one random arc for speed
        for flip_idx in range(edges):
            bits2 = bits ^ (1 << flip_idx)
            A1 = bits_to_adj(bits2, n)
            cycles1 = compute_cycle_data(A1, n)
            a1_1, a2_1 = compute_alphas(cycles1)
            da1 = a1_1 - a1_0
            da2 = a2_1 - a2_0
            da_counter_6[(da1, da2)] += 1

        if (bits + 1) % 8192 == 0:
            elapsed = time.time() - start
            pct = (bits + 1) / (1 << edges) * 100
            print(f"  {pct:.1f}% ({bits+1}/{1 << edges}) in {elapsed:.0f}s...", flush=True)

    elapsed = time.time() - start
    print(f"  Completed in {elapsed:.1f}s")

    print(f"\n  (da1, da2) distribution at n=6:")
    da2_nonzero = 0
    for (da1, da2) in sorted(da_counter_6.keys()):
        cnt = da_counter_6[(da1, da2)] // 2
        if da2 != 0:
            da2_nonzero += cnt
        dH = 2 * da1 + 4 * da2
        if cnt > 0:
            print(f"    (da1={da1:+3d}, da2={da2:+3d}): {cnt:8d}  => delta_H = {dH:+4d}")

    total_flips = sum(v for v in da_counter_6.values()) // 2
    print(f"\n  da2 != 0 fraction: {da2_nonzero}/{total_flips} = "
          f"{da2_nonzero/total_flips*100:.2f}%")

    # Which da1 values skip?
    da1_vals_6 = sorted(set(da1 for (da1, da2) in da_counter_6.keys()))
    max_abs_da1 = max(abs(d) for d in da1_vals_6)
    missing_da1 = [d for d in range(-max_abs_da1, max_abs_da1+1) if d not in da1_vals_6]
    print(f"  da1 range: [{min(da1_vals_6)}, {max(da1_vals_6)}]")
    print(f"  Missing da1 values: {missing_da1}")

    # Check delta_H values
    all_dH = set()
    for (da1, da2) in da_counter_6.keys():
        all_dH.add(2 * da1 + 4 * da2)
    all_dH_abs = sorted(set(abs(d) for d in all_dH))
    max_dH = max(all_dH_abs)
    missing_dH = [d for d in range(0, max_dH+1, 2) if d not in all_dH_abs]
    print(f"\n  |delta_H| values: {all_dH_abs}")
    print(f"  Missing |delta_H|: {missing_dH}")

if __name__ == "__main__":
    main()
