"""
h_landscape_deep.py -- kind-pasteur-2026-03-14-S71
Deep investigation of the H-landscape on the tournament hypercube.

KEY DISCOVERY FROM PREVIOUS SCRIPT:
  At n=4,5: the H-landscape has a UNIQUE local maximum (the global max)
  This means gradient ascent ALWAYS finds the H-maximizer!

QUESTIONS:
1. Does this hold at n=6,7? (exhaustive n=6, sampled n=7)
2. What is the "basin of attraction" of the global maximum?
3. How many gradient ascent steps needed from random start?
4. Is the landscape "unimodal" (no local maxima other than global)?
5. Connection to Condorcet/social choice: monotonicity of H under flips?
6. Does the H-landscape have any FLAT regions (plateaus)?
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys, random

sys.stdout.reconfigure(encoding='utf-8')

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def gradient_ascent(bits, n, m):
    """
    Greedy gradient ascent on H-landscape.
    At each step, flip the arc that maximizes H.
    Return: (final_bits, final_H, steps, path_of_H_values)
    """
    A = bits_to_adj(bits, n)
    H = compute_H_dp(A, n)
    path = [H]
    steps = 0

    while True:
        best_nbr_bits = None
        best_nbr_H = H

        for arc_bit in range(m):
            nbr_bits = bits ^ (1 << arc_bit)
            A_nbr = bits_to_adj(nbr_bits, n)
            H_nbr = compute_H_dp(A_nbr, n)
            if H_nbr > best_nbr_H:
                best_nbr_H = H_nbr
                best_nbr_bits = nbr_bits

        if best_nbr_bits is None:
            # No improvement — at local maximum
            break

        bits = best_nbr_bits
        H = best_nbr_H
        A = bits_to_adj(bits, n)
        path.append(H)
        steps += 1

    return bits, H, steps, path

def main():
    print("=" * 70)
    print("H-LANDSCAPE DEEP INVESTIGATION")
    print("kind-pasteur-2026-03-14-S71")
    print("=" * 70)

    # ========================================
    # PART 1: Exhaustive at n=6
    # ========================================
    print("\n" + "=" * 70)
    print("PART 1: EXHAUSTIVE LOCAL EXTREMA AT n=6")
    print("  n=6: 2^15 = 32768 tournaments, 15-dim hypercube")
    print("=" * 70)

    n = 6
    m = n * (n-1) // 2
    total = 2**m

    local_max_count = 0
    local_min_count = 0
    max_H_vals = set()
    min_H_vals = set()
    plateau_count = 0  # tournaments where some neighbor has same H

    for bits in range(total):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)

        is_max = True
        is_min = True
        has_flat = False

        for arc_bit in range(m):
            nbr_bits = bits ^ (1 << arc_bit)
            A_nbr = bits_to_adj(nbr_bits, n)
            H_nbr = compute_H_dp(A_nbr, n)
            if H_nbr > H:
                is_max = False
            if H_nbr < H:
                is_min = False
            if H_nbr == H:
                has_flat = True

        if is_max:
            local_max_count += 1
            max_H_vals.add(H)
        if is_min:
            local_min_count += 1
            min_H_vals.add(H)
        if has_flat:
            plateau_count += 1

    print(f"  Total: {total}")
    print(f"  Local maxima: {local_max_count} ({100*local_max_count/total:.2f}%)")
    print(f"  Local minima: {local_min_count} ({100*local_min_count/total:.2f}%)")
    print(f"  Plateau (some flat neighbor): {plateau_count} ({100*plateau_count/total:.2f}%)")
    print(f"  H at local maxima: {sorted(max_H_vals)}")
    print(f"  H at local minima: {sorted(min_H_vals)}")
    print(f"  UNIQUE local max? {len(max_H_vals) == 1}")

    # ========================================
    # PART 2: Gradient ascent statistics
    # ========================================
    print("\n" + "=" * 70)
    print("PART 2: GRADIENT ASCENT FROM RANDOM STARTS")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n--- n = {n} ---")
        m = n * (n-1) // 2
        total = 2**m
        random.seed(42)

        num_trials = min(500, total)
        steps_list = []
        final_H_list = []

        for trial in range(num_trials):
            if num_trials < total:
                bits = random.randint(0, total - 1)
            else:
                bits = trial
            _, final_H, steps, path = gradient_ascent(bits, n, m)
            steps_list.append(steps)
            final_H_list.append(final_H)

        print(f"  {num_trials} gradient ascent runs:")
        print(f"  Final H distribution: {dict(Counter(final_H_list).most_common(5))}")
        print(f"  All reach global max? {len(set(final_H_list)) == 1}")
        print(f"  Steps: mean={np.mean(steps_list):.2f}, max={max(steps_list)}, "
              f"median={np.median(steps_list):.0f}")
        print(f"  Steps distribution: {dict(Counter(steps_list).most_common(10))}")

    # ========================================
    # PART 3: Sampled at n=7
    # ========================================
    print("\n" + "=" * 70)
    print("PART 3: SAMPLED LOCAL EXTREMA AT n=7")
    print("  Check: does gradient ascent always reach H=189?")
    print("=" * 70)

    n = 7
    m = n * (n-1) // 2
    random.seed(123)

    num_trials = 200
    final_H_list = []
    steps_list = []
    local_max_H = set()

    for trial in range(num_trials):
        bits = random.randint(0, 2**m - 1)
        final_bits, final_H, steps, path = gradient_ascent(bits, n, m)
        final_H_list.append(final_H)
        steps_list.append(steps)
        local_max_H.add(final_H)

    print(f"  {num_trials} gradient ascent runs:")
    print(f"  Final H values: {sorted(local_max_H)}")
    print(f"  Final H distribution: {dict(Counter(final_H_list).most_common(10))}")
    print(f"  All reach H=189? {set(final_H_list) == {189}}")
    print(f"  Steps: mean={np.mean(steps_list):.2f}, max={max(steps_list)}")

    # ========================================
    # PART 4: Diameter and connectivity
    # ========================================
    print("\n" + "=" * 70)
    print("PART 4: H-FIBER CONNECTIVITY")
    print("  Are all tournaments with H=k connected by arc flips")
    print("  that preserve or increase H?")
    print("=" * 70)

    n = 5
    m = n * (n-1) // 2
    total = 2**m

    # BFS from any maximizer to find how many tournaments are "reachable"
    # via non-decreasing H paths
    print(f"\n--- n = {n} ---")

    # Find all maximizers
    max_H = 0
    for bits in range(total):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        max_H = max(max_H, H)

    maximizers = set()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        if H == max_H:
            maximizers.add(bits)

    print(f"  Global max H = {max_H}, # maximizers = {len(maximizers)}")

    # From each H value, how many steps DOWN (non-increasing) to reach H=1?
    # And from H=1, how many steps UP to reach H=max?
    H_map = {}
    for bits in range(total):
        A = bits_to_adj(bits, n)
        H_map[bits] = compute_H_dp(A, n)

    # Count H-monotone paths
    print(f"\n  Gradient descent from maximizer to H=1:")
    for start_bits in list(maximizers)[:3]:
        # Random walk always choosing decreasing H
        bits = start_bits
        H = H_map[bits]
        path = [H]
        while H > 1:
            candidates = []
            for arc_bit in range(m):
                nbr = bits ^ (1 << arc_bit)
                if H_map[nbr] < H:
                    candidates.append((H_map[nbr], nbr))
            if not candidates:
                break
            # Choose smallest H neighbor
            candidates.sort()
            H = candidates[0][0]
            bits = candidates[0][1]
            path.append(H)
        print(f"    Path: {path}")

    # ========================================
    # PART 5: Score sequence and gradient
    # ========================================
    print("\n" + "=" * 70)
    print("PART 5: ARC-FLIP GRADIENT vs SCORE REGULARIZATION")
    print("  Does maximizing H = regularizing scores?")
    print("=" * 70)

    n = 5
    m = n * (n-1) // 2
    random.seed(42)

    for trial in range(5):
        bits = random.randint(0, 2**m - 1)
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        scores = sorted([sum(A[i]) for i in range(n)])
        score_var = np.var(scores)

        print(f"\n  Trial {trial}: initial H={H}, scores={scores}, var={score_var:.2f}")

        _, final_H, steps, path = gradient_ascent(bits, n, m)
        A_final = bits_to_adj(_, n)
        final_scores = sorted([sum(A_final[i]) for i in range(n)])
        final_var = np.var(final_scores)

        print(f"    Final H={final_H}, scores={final_scores}, var={final_var:.2f}, steps={steps}")
        print(f"    Score variance: {score_var:.2f} -> {final_var:.2f} (change={final_var-score_var:+.2f})")

    # ========================================
    # PART 6: The "H=7 barrier" in the landscape
    # ========================================
    print("\n" + "=" * 70)
    print("PART 6: H=7 GAP IN THE LANDSCAPE")
    print("  H=7 is impossible. How does the landscape 'skip' it?")
    print("  At n=5: H goes 1,3,5 then jumps to 9,11,13,15")
    print("  Are all H=5 tournaments adjacent to H=9 tournaments?")
    print("=" * 70)

    n = 5
    m = n * (n-1) // 2

    # For each H=5 tournament, what are the neighbor H values?
    neighbors_of_5 = Counter()
    count_5 = 0
    for bits in range(2**m):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        if H == 5:
            count_5 += 1
            for arc_bit in range(m):
                nbr = bits ^ (1 << arc_bit)
                A_nbr = bits_to_adj(nbr, n)
                H_nbr = compute_H_dp(A_nbr, n)
                neighbors_of_5[H_nbr] += 1

    print(f"\n  n=5: {count_5} tournaments with H=5")
    print(f"  Neighbor H values from H=5:")
    for H_nbr in sorted(neighbors_of_5.keys()):
        avg = neighbors_of_5[H_nbr] / count_5
        print(f"    H={H_nbr}: {neighbors_of_5[H_nbr]} total ({avg:.1f} per H=5 tournament)")

    # Same for H=9
    neighbors_of_9 = Counter()
    count_9 = 0
    for bits in range(2**m):
        A = bits_to_adj(bits, n)
        H = compute_H_dp(A, n)
        if H == 9:
            count_9 += 1
            for arc_bit in range(m):
                nbr = bits ^ (1 << arc_bit)
                A_nbr = bits_to_adj(nbr, n)
                H_nbr = compute_H_dp(A_nbr, n)
                neighbors_of_9[H_nbr] += 1

    print(f"\n  {count_9} tournaments with H=9")
    print(f"  Neighbor H values from H=9:")
    for H_nbr in sorted(neighbors_of_9.keys()):
        avg = neighbors_of_9[H_nbr] / count_9
        print(f"    H={H_nbr}: {neighbors_of_9[H_nbr]} total ({avg:.1f} per H=9 tournament)")

    # KEY: H=5 -> H=9 requires jumping by +4. Is this always via a single flip?
    direct_5_to_9 = 0
    for bits in range(2**m):
        if H_map.get(bits, compute_H_dp(bits_to_adj(bits, n), n)) == 5:
            for arc_bit in range(m):
                nbr = bits ^ (1 << arc_bit)
                if compute_H_dp(bits_to_adj(nbr, n), n) == 9:
                    direct_5_to_9 += 1
    print(f"\n  Direct H=5 -> H=9 single-flip edges: {direct_5_to_9}")
    print(f"  H=7 truly SKIPPED in the landscape (no tournament has H=7)")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    main()
