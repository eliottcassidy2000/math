#!/usr/bin/env python3
"""
fiber_fractions_s20aw.py -- kind-pasteur-2026-03-22-S20aw

FIBER FRACTIONS: The exact arithmetic of score fibers.

The score map pi: {0,1}^m -> score sequences is a FIBER BUNDLE.
Each fiber pi^{-1}(s) contains all tournaments with score sequence s.

Key quantities to compute exactly:
1. FIBER SIZE: |pi^{-1}(s)| = number of tournaments with score s
2. FIBER FRACTION: |pi^{-1}(s)| / 2^m = probability of score s
3. WITHIN-FIBER EDGE FRACTION: what fraction of arc flips stay in the fiber
4. FIBER CONNECTIVITY: how many components within each fiber
5. H VARIATION WITHIN FIBER: range and variance of H within each fiber
6. FIBER MULTIPLICITY: how many iso classes per fiber

Building on opus S191's findings:
- Within-fiber fractions: 1/2, 3/8, 5/16 at n=3,4,5
- Score-preserving flips need |s_i - s_j| <= 1
- Disconnected fibers are the NORM (only 3/9 connected at n=5)
- The (0,2,2,2,4) fiber has 0 within-fiber edges (totally disconnected)

Author: kind-pasteur-2026-03-22-S20aw
"""
import sys
import numpy as np
from math import comb, factorial, log2
from collections import defaultdict
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

print("=" * 70)
print("  FIBER FRACTIONS: EXACT ARITHMETIC OF SCORE FIBERS")
print("=" * 70)

# ================================================================
# EXACT COMPUTATION AT n=3,4,5,6
# ================================================================
for n in [3, 4, 5, 6]:
    print(f"\n{'='*70}")
    print(f"  n = {n}")
    print(f"{'='*70}\n")

    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    total = 2**m

    # Enumerate all tournaments, group by score
    fibers = defaultdict(list)  # score -> list of (bits, H)
    for bits in range(total):
        A = np.zeros((n,n), dtype=np.int8)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        score = tuple(sorted(A.sum(axis=1).astype(int)))
        H = count_hp(A, n)
        fibers[score].append((bits, H))

    # Count within-fiber edges
    total_flips = 0
    within_flips = 0
    fiber_data = []

    bits_to_score = {}
    for score, members in fibers.items():
        for bits, H in members:
            bits_to_score[bits] = score

    for bits in range(total):
        for k in range(m):
            total_flips += 1
            nb = bits ^ (1 << k)
            if bits_to_score[bits] == bits_to_score[nb]:
                within_flips += 1

    within_frac = Fraction(within_flips, total_flips)

    print(f"  Total tournaments: {total}")
    print(f"  Score sequences: {len(fibers)}")
    print(f"  Within-fiber flip fraction: {within_flips}/{total_flips} = {within_frac} = {float(within_frac):.6f}")
    print()

    # Detailed fiber table
    print(f"  {'Score':>25s} {'Size':>6s} {'Frac':>10s} {'#H-vals':>8s} {'H range':>12s} {'Mean H':>8s} {'Within%':>8s}")

    for score in sorted(fibers.keys()):
        members = fibers[score]
        size = len(members)
        frac = Fraction(size, total)
        H_vals = sorted(set(h for _, h in members))
        mean_H = sum(h for _, h in members) / size

        # Within-fiber edges for this specific fiber
        member_set = set(b for b, _ in members)
        fiber_within = 0
        fiber_total = 0
        for bits, _ in members:
            for k in range(m):
                fiber_total += 1
                nb = bits ^ (1 << k)
                if nb in member_set:
                    fiber_within += 1

        within_pct = 100 * fiber_within / fiber_total if fiber_total > 0 else 0

        H_range_str = f"[{min(H_vals)},{max(H_vals)}]" if len(H_vals) > 1 else f"{H_vals[0]}"

        fiber_data.append({
            'score': score, 'size': size, 'frac': frac,
            'H_vals': H_vals, 'mean_H': mean_H, 'within_pct': within_pct
        })

        print(f"  {str(list(score)):>25s} {size:>6d} {str(frac):>10s} {len(H_vals):>8d} {H_range_str:>12s} {mean_H:>8.1f} {within_pct:>7.1f}%")

    # ================================================================
    # THE FIBER FRACTION FORMULA
    # ================================================================
    print(f"\n  FIBER SIZE FORMULA:")

    # The number of tournaments with score sequence s = (s_0, ..., s_{n-1})
    # is a multinomial-type count. For sorted score s, the count is:
    # |fiber(s)| = n! / |Aut(s)| * T(s)
    # where T(s) is the number of tournaments on {1,...,n} with scores s
    # (up to the labeling determined by s).
    # Actually: |fiber(s)| = (multinomial coefficient) * (tournament count).
    # The multinomial part accounts for how many ways to assign the scores to vertices.

    # Simpler: the multiplicity of score s in the sorted sense is:
    # |fiber(s)| = n! / prod(k_i!) * t(s)
    # where k_i is the multiplicity of score value i, and t(s) is the
    # number of labeled tournaments with THAT specific vertex-to-score assignment.

    for fd in fiber_data:
        score = fd['score']
        # Compute the score multiplicity: how many ways to assign these scores to vertices
        from collections import Counter
        multiplicities = Counter(score)
        multinomial = factorial(n)
        for cnt in multiplicities.values():
            multinomial //= factorial(cnt)
        # The number of tournaments per assignment = size / multinomial
        per_assignment = Fraction(fd['size'], multinomial)
        print(f"    {list(score)}: size={fd['size']}, mult_coeff={multinomial}, per_assign={per_assignment}")

    # ================================================================
    # THE S2 -> FIBER SIZE MAPPING
    # ================================================================
    print(f"\n  S2 -> FIBER SIZE:")
    s2_to_sizes = defaultdict(list)
    for fd in fiber_data:
        s2 = sum(s*s for s in fd['score']) - comb(n, 2)
        s2_to_sizes[s2].append((fd['score'], fd['size']))

    for s2 in sorted(s2_to_sizes.keys()):
        entries = s2_to_sizes[s2]
        total_size = sum(sz for _, sz in entries)
        print(f"    S2-C(n,2)={s2}: {len(entries)} score seq(s), total size={total_size}")

    # ================================================================
    # CROSS-FIBER TRANSITION MATRIX
    # ================================================================
    print(f"\n  CROSS-FIBER TRANSITION MATRIX:")
    print(f"  (Fraction of flips from fiber s that land in fiber s')")

    score_list = sorted(fibers.keys())
    n_fibers = len(score_list)
    score_idx = {s: i for i, s in enumerate(score_list)}

    # Build transition counts
    T = np.zeros((n_fibers, n_fibers), dtype=int)
    for score, members in fibers.items():
        i = score_idx[score]
        member_set = set(b for b, _ in members)
        for bits, _ in members:
            for k in range(m):
                nb = bits ^ (1 << k)
                nb_score = bits_to_score[nb]
                j = score_idx[nb_score]
                T[i][j] += 1

    # Normalize by row
    print(f"\n  Transition fractions (rows = source fiber, cols = dest fiber):")
    print(f"  {'src\\dst':>10s}", end="")
    for j in range(n_fibers):
        print(f" {str(list(score_list[j])):>12s}", end="")
    print()

    for i in range(n_fibers):
        row_sum = T[i].sum()
        print(f"  {str(list(score_list[i])):>10s}", end="")
        for j in range(n_fibers):
            if T[i][j] > 0:
                frac = Fraction(int(T[i][j]), int(row_sum))
                print(f" {str(frac):>12s}", end="")
            else:
                print(f" {'0':>12s}", end="")
        print()

    # ================================================================
    # THE EXACT WITHIN-FIBER FRACTION
    # ================================================================
    print(f"\n  WITHIN-FIBER FRACTION: {within_frac}")
    if n <= 5:
        # Check if it matches a simple formula
        # From opus S191: 1/2, 3/8, 5/16
        # Numerators: 1, 3, 5 (odd numbers)
        # Denominators: 2, 8, 16 = 2^1, 2^3, 2^4
        # Actually: 2, 8, 16 = 2, 2^3, 2^4. Hmm.
        # Let me check: 1/2, 3/8, 5/16
        # = 8/16, 6/16, 5/16? No that's wrong.
        # 1/2 = 4/8, 3/8, ...
        # Denominators: 2 = 2^1, 8 = 2^3, 16 = 2^4
        # These are 2^{C(n-1,2)}: C(2,2)=1, C(3,2)=3, C(4,2)=6
        # 2^1, 2^3, 2^6 = 2, 8, 64. But denominator at n=5 is 16, not 64.
        # So the 2^{C(n-1,2)} pattern from opus S181 was for self-loops, not within-fiber.
        pass

print(f"\n{'='*70}")
print(f"  WITHIN-FIBER FRACTION SEQUENCE")
print(f"{'='*70}\n")

# Print the sequence we computed
print(f"  n=3: within-fiber fraction = 1/2")
print(f"  n=4: within-fiber fraction = 3/8")
print(f"  n=5: within-fiber fraction = 5/16")
# n=6 computed above
