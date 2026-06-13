#!/usr/bin/env python3
"""
h63_gap_check.py -- Check whether H=63 is a permanent gap in the H-spectrum.

Context:
  H(T) = number of directed Hamiltonian paths in tournament T.
  By OCF: H(T) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

  Known gaps: H=7 (proved impossible for all n), H=21 (strongly conjectured).
  At n=7: H=63 is absent (gap between 61 and 65 in the exhaustive spectrum).
  Pattern: 7, 21=7*3, 63=7*9. Also 7=2^3-1, 63=2^6-1.

Approach:
  1. OCF decomposition analysis for H=63.
  2. Sample random n=8 tournaments (500k) and n=9 tournaments (200k).
  3. Also do near-transitive perturbations.
  4. Track H=21 as a control (should remain absent).
  5. Report distribution of H values near 63.
  6. Analyze the 2^k-1 pattern among gaps.

The exhaustive n=7 spectrum is already verified (77 distinct values,
H=63 absent -- see OPEN-QUESTIONS.md).
"""

import random
import time
import sys
from collections import Counter

# ============================================================
# Core: Hamiltonian path count via bitmask DP (optimized with dict)
# ============================================================

def count_H_dp(adj, n):
    """Count Hamiltonian paths via bitmask DP using sparse dict."""
    full = (1 << n) - 1
    # Use dict for sparse storage: dp[(mask, v)] = count
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c

    return sum(dp.get((full, v), 0) for v in range(n))


def count_H_dp_fast(adj_rows, n):
    """Faster DP using precomputed adjacency as bitmasks per row."""
    full = (1 << n) - 1
    # dp[mask] is a list of length n
    # Use flat array approach
    size = 1 << n
    # dp[mask * n + v] = count
    dp = [0] * (size * n)

    for v in range(n):
        dp[(1 << v) * n + v] = 1

    for mask in range(1, size):
        base = mask * n
        for v in range(n):
            c = dp[base + v]
            if c == 0:
                continue
            # neighbors of v not in mask
            targets = adj_rows[v] & ~mask & full
            u = 0
            t = targets
            while t:
                if t & 1:
                    dp[(mask | (1 << u)) * n + u] += c
                t >>= 1
                u += 1

    return sum(dp[full * n + v] for v in range(n))


def random_tournament_bitmask(n):
    """Generate random tournament, return adjacency as bitmask rows."""
    adj_rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                adj_rows[i] |= (1 << j)
            else:
                adj_rows[j] |= (1 << i)
    return adj_rows


def near_transitive_bitmask(n, num_flips):
    """Transitive tournament with random edge flips, as bitmask rows."""
    adj_rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj_rows[i] |= (1 << j)

    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    to_flip = random.sample(edges, min(num_flips, len(edges)))
    for i, j in to_flip:
        if adj_rows[i] & (1 << j):
            adj_rows[i] &= ~(1 << j)
            adj_rows[j] |= (1 << i)
        else:
            adj_rows[j] &= ~(1 << i)
            adj_rows[i] |= (1 << j)
    return adj_rows


# ============================================================
# Part 1: OCF decomposition analysis
# ============================================================

def analyze_ocf_decompositions():
    print("=" * 70)
    print("OCF DECOMPOSITION ANALYSIS FOR H=63")
    print("=" * 70)
    print()
    print("H(T) = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4 + 32*a5 + ...")
    print("H = 63 requires: a1 + 2*a2 + 4*a3 + 8*a4 + 16*a5 = 31")
    print()

    solutions = []
    for a5 in range(2):
        rem4 = 31 - 16 * a5
        for a4 in range(rem4 // 8 + 1):
            rem3 = rem4 - 8 * a4
            for a3 in range(rem3 // 4 + 1):
                rem2 = rem3 - 4 * a3
                for a2 in range(rem2 // 2 + 1):
                    a1 = rem2 - 2 * a2
                    solutions.append((a1, a2, a3, a4, a5))

    print(f"Number of non-negative integer solutions: {len(solutions)}")
    print()
    # Just show a few representative ones
    print("Sample solutions (a1, a2, a3, a4, a5):")
    for sol in solutions[:10]:
        print(f"  {sol}")
    if len(solutions) > 10:
        print(f"  ... ({len(solutions) - 10} more)")
    print()

    # Key insight: binary representation of 31 = 11111 in binary
    print("Key: 31 = 11111 in binary, so the 'all ones' solution is")
    print("  (1, 1, 1, 1, 1) -> sum = 1+2+4+8+16 = 31  [valid]")
    print()
    print("Note: 63 = 2^6 - 1 = 111111 in binary")
    print("      H=63 = 1 + 62 = 1 + 2*31")
    print("      The '1' comes from the transitive tournament's single path")
    print()


# ============================================================
# Part 2: Sampling
# ============================================================

def sample_tournaments(n, num_random, num_near_trans):
    print(f"\n{'=' * 70}")
    print(f"SAMPLING n={n}: {num_random} random + {num_near_trans} near-transitive")
    print(f"{'=' * 70}")

    h_counter = Counter()
    t0 = time.time()

    # Random sampling
    for trial in range(num_random):
        adj = random_tournament_bitmask(n)
        h = count_H_dp_fast(adj, n)
        h_counter[h] += 1

        if (trial + 1) % 25000 == 0:
            elapsed = time.time() - t0
            rate = (trial + 1) / elapsed
            print(f"  Random: {trial+1}/{num_random}, {rate:.0f}/sec", flush=True)

    t1 = time.time()
    print(f"  Random phase: {t1-t0:.1f}s")

    # Near-transitive sampling (targets low H values)
    for trial in range(num_near_trans):
        num_flips = random.randint(1, n * (n - 1) // 4)
        adj = near_transitive_bitmask(n, num_flips)
        h = count_H_dp_fast(adj, n)
        h_counter[h] += 1

        if (trial + 1) % 25000 == 0:
            elapsed = time.time() - t1
            rate = (trial + 1) / elapsed
            print(f"  Near-trans: {trial+1}/{num_near_trans}, {rate:.0f}/sec",
                  flush=True)

    elapsed = time.time() - t0
    total = num_random + num_near_trans
    print(f"  Total: {total} samples in {elapsed:.1f}s ({total/elapsed:.0f}/sec)")

    # Results
    print()
    for target in [7, 21, 63]:
        if target in h_counter:
            print(f"  H={target}: *** FOUND *** Count = {h_counter[target]}")
        else:
            print(f"  H={target}: not found")

    print()
    print("  Distribution near H=63:")
    for h in range(53, 77, 2):
        count = h_counter.get(h, 0)
        bar = "#" * min(count, 50)
        print(f"    H={h:3d}: {count:6d}  {bar}")

    all_h = sorted(h_counter.keys())
    print(f"\n  Distinct H values: {len(all_h)}, range [{min(all_h)}, {max(all_h)}]")

    # Show all gaps in [1, 100] that are odd and not found
    small_gaps = [h for h in range(1, min(102, max(all_h)), 2) if h not in h_counter]
    if small_gaps:
        print(f"  Odd gaps in [1,101]: {small_gaps}")

    return h_counter


# ============================================================
# Part 3: Pattern analysis
# ============================================================

def analyze_pattern():
    # Known n=7 gaps from exhaustive enumeration
    n7_spectrum = [1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33,
                   35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61,
                   65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91,
                   93, 95, 97, 99, 101, 103, 105, 109, 111, 113, 115, 117,
                   121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141,
                   143, 145, 147, 151, 153, 155, 157, 159, 171, 175, 189]

    n7_set = set(n7_spectrum)
    max_h = max(n7_spectrum)

    # All odd gaps
    gaps = sorted(h for h in range(1, max_h + 1, 2) if h not in n7_set)

    print()
    print("=" * 70)
    print("PATTERN ANALYSIS OF n=7 GAPS")
    print("=" * 70)
    print()
    print(f"All {len(gaps)} gaps at n=7 (odd values in [1,{max_h}] not achieved):")
    print(f"  {gaps}")
    print()

    # Categorize gaps
    print("Gap categorization:")
    for g in gaps:
        is_mersenne = (g & (g + 1)) == 0  # 2^k - 1
        div7 = g % 7 == 0

        # Factor
        factors = []
        temp = g
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
            while temp % p == 0:
                factors.append(p)
                temp //= p
        if temp > 1:
            factors.append(temp)

        tag = ""
        if is_mersenne:
            k = g.bit_length()
            tag += f" [= 2^{k}-1, Mersenne]"
        if div7:
            tag += f" [= 7*{g//7}]"

        print(f"  {g:3d} = {'*'.join(map(str,factors)):15s}{tag}")

    print()
    print("Multiples of 7 (odd) in [1,189]:")
    for m in range(7, 190, 14):  # 7, 21, 35, ...
        status = "GAP" if m not in n7_set else "achieved"
        is_mersenne = (m & (m + 1)) == 0
        extra = " [Mersenne]" if is_mersenne else ""
        print(f"  7*{m//7:2d} = {m:3d}: {status:10s}{extra}")

    print()
    print("Mersenne numbers 2^k-1 in range:")
    for k in range(1, 9):
        m = (1 << k) - 1
        if m <= max_h:
            status = "GAP" if m not in n7_set else "achieved"
        else:
            status = "out of range"
        print(f"  2^{k}-1 = {m:3d}: {status}")

    print()
    print("STRUCTURAL OBSERVATIONS:")
    print()
    print("1. H=7 = 2^3-1: proved permanent gap (no OCF decomposition works)")
    print("2. H=21 = 3*7: strongly conjectured permanent gap")
    print("3. H=63 = 9*7 = 2^6-1: absent at n=7, status unknown")
    print()
    print("The '7 family' hypothesis: gaps at 7, 21=3*7, 63=9*7, ...")
    print("  But 35=5*7 is ACHIEVED at n=7, so not all multiples of 7 are gaps!")
    print("  Also 49=7*7 is achieved, 77=11*7 is achieved, etc.")
    print()
    print("The 'Mersenne' hypothesis: gaps at 2^k-1")
    print("  7 = 2^3-1: GAP")
    print("  63 = 2^6-1: GAP (at n=7)")
    print("  But 1=2^1-1: achieved, 3=2^2-1: achieved, 15=2^4-1: achieved,")
    print("  31=2^5-1: achieved. So not all Mersenne numbers are gaps.")
    print()
    print("The 'divisible by 7 AND Mersenne-like' pattern:")
    print("  7, 21, 63 are gaps. They share: all are divisible by 7.")
    print("  21 = 3*(2^3-1), 63 = 9*(2^3-1) = (2^6-1)")
    print("  BUT: 35=5*7 and 49=7*7 are NOT gaps.")
    print("  Perhaps the constraint is more nuanced (OCF binary structure).")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    random.seed(42)

    # Part 1: OCF analysis
    analyze_ocf_decompositions()

    # Part 2: Pattern analysis (no computation needed)
    analyze_pattern()

    # Part 3: Sampling at n=8
    h8 = sample_tournaments(8, 500000, 100000)

    # Part 4: Sampling at n=9 (slower, fewer samples)
    h9 = sample_tournaments(9, 200000, 50000)

    # Part 5: Summary
    print()
    print("=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print()

    found_63_n8 = 63 in h8
    found_63_n9 = 63 in h9
    found_21_n8 = 21 in h8
    found_21_n9 = 21 in h9

    print("H=7:  permanent gap (proved, THM-029)")

    if found_21_n8 or found_21_n9:
        print(f"H=21: *** FOUND at n={'8' if found_21_n8 else '9'} ***")
    else:
        print(f"H=21: absent in all samples (control: consistent with conjecture)")

    if found_63_n8 or found_63_n9:
        n_found = '8' if found_63_n8 else '9'
        c = h8.get(63, 0) + h9.get(63, 0)
        print(f"H=63: *** FOUND at n={n_found}, count={c} *** NOT a permanent gap")
    else:
        print(f"H=63: absent in 850k samples at n=8,9")
        print(f"       Strong evidence for permanent gap (like H=21)")

    # Check what other n=7 gaps persist
    n7_gaps = {7, 21, 63, 107, 119, 149,
               161, 163, 165, 167, 169, 173,
               177, 179, 181, 183, 185, 187}
    print()
    print("Status of n=7 gaps at n=8,9 (from sampling):")
    for g in sorted(n7_gaps):
        in8 = g in h8
        in9 = g in h9
        if in8 or in9:
            status = f"FOUND at n={'8' if in8 else '9'}"
        else:
            status = "still absent"
        print(f"  H={g:3d}: {status}")
