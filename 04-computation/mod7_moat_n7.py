"""
mod7_moat_n7.py -- kind-pasteur-2026-03-14-S67

At n=6, H mod 7 = 0 is COMPLETELY impossible (0/32768 tournaments).
Root cause: alpha_1 + 2*alpha_2 never equals 3 mod 7.

At n=7, H=189 = 27*7 IS achievable. So mod 7 barrier breaks.
Questions:
1. Which specific multiples of 7 are achievable at n=7?
2. What is the (alpha_1, alpha_2) mod 7 distribution at n=7?
3. Is there a refined modular constraint blocking 7, 21, 63 specifically?

Note: at n=7, H = 1 + 2*a1 + 4*a2 + 8*a3 (alpha_3 can be nonzero).
Max IS size = 2 at n=7 (two disjoint 3-cycles = 6 verts, or 3+5=8>7, or 3+3=6<=7).
Wait: at n=7, can we have three disjoint 3-cycles? Need 9 > 7. No.
Can we have (3-cycle, 3-cycle)? Need 6 <= 7. Yes.
But 5-cycle + 3-cycle = 8 > 7? Wait, they need disjoint vertex sets.
5+3 = 8 > 7. Can't fit. But two 3-cycles: 3+3=6 <= 7. CAN fit.
So max IS size at n=7 = 2, and alpha_3 = 0. H = 1 + 2*a1 + 4*a2.
But wait, I need to count WEIGHTED independent sets, where each "vertex"
in Omega is a directed cycle weighted by its count.
Actually alpha_k counts products of weights for each independent k-set.
For IS of size 2: pairs of disjoint cycles, weighted by product of their counts.
So alpha_2 is nonzero at n=7, alpha_3=0 (can't fit 3 pairwise disjoint odd cycles in 7 vertices).
Wait - could we have three disjoint 3-cycles in 7 vertices? No: 3*3=9>7.
Two disjoint 3-cycles in 7 vertices? 2*3=6<=7. Yes.
So alpha_2 counts (among other things) pairs of disjoint 3-cycles.
alpha_3 = 0 at n=7.

H = 1 + 2*alpha_1 + 4*alpha_2 at n=7.
"""

import numpy as np
from itertools import combinations
from collections import Counter
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

def compute_alpha_n7(A, n):
    """Compute alpha_1, alpha_2 for n=7 (alpha_3=0)."""
    cycles = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt, size))

    alpha_1 = sum(cnt for _, cnt, _ in cycles)

    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) == 0:
                alpha_2 += cycles[i][1] * cycles[j][1]

    return alpha_1, alpha_2, cycles

def main():
    print("=" * 70)
    print("MOD 7 MOAT ANALYSIS AT n=7 (SAMPLING)")
    print("=" * 70)

    n = 7
    rng = np.random.default_rng(2026_03_14_97)

    H_counter = Counter()
    pair_counter = Counter()  # (a1 mod 7, a2 mod 7)
    sum_mod7 = Counter()  # (a1 + 2*a2) mod 7
    multiples_of_7 = Counter()

    # Track (a1, a2) pairs that give H divisible by 7
    h_mod7_examples = {}

    num_samples = 300000
    start = time.time()

    for trial in range(num_samples):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        a1, a2, cycles = compute_alpha_n7(A, n)
        H = 1 + 2*a1 + 4*a2
        H_counter[H] += 1
        pair_counter[(a1 % 7, a2 % 7)] += 1
        s = (a1 + 2*a2) % 7
        sum_mod7[s] += 1

        if H % 7 == 0:
            multiples_of_7[H] += 1
            if H not in h_mod7_examples:
                h_mod7_examples[H] = (a1, a2)

        if (trial+1) % 50000 == 0:
            elapsed = time.time() - start
            print(f"  {trial+1}/{num_samples} ({(trial+1)/elapsed:.0f}/sec)...", flush=True)

    elapsed = time.time() - start
    print(f"\nCompleted {num_samples} samples in {elapsed:.1f}s")
    print(f"Distinct H values: {len(H_counter)}")

    # 1. H mod 7 distribution
    print(f"\n{'='*70}")
    print("H mod 7 DISTRIBUTION AT n=7")
    print("=" * 70)
    hmod7 = Counter()
    for h, cnt in H_counter.items():
        hmod7[h % 7] += cnt
    for r in range(7):
        pct = 100 * hmod7[r] / num_samples
        print(f"  H mod 7 = {r}: {hmod7[r]:7d} ({pct:.2f}%)")

    # 2. Which multiples of 7 are achievable?
    print(f"\n{'='*70}")
    print("ACHIEVABLE MULTIPLES OF 7 AT n=7")
    print("=" * 70)
    if multiples_of_7:
        for h in sorted(multiples_of_7.keys()):
            ex = h_mod7_examples.get(h, ("?","?"))
            print(f"  H={h:4d} = 7*{h//7:3d}: {multiples_of_7[h]:6d} hits  (e.g. a1={ex[0]}, a2={ex[1]})")
    else:
        print("  NONE found!")

    # Missing multiples of 7
    max_H = max(H_counter.keys())
    print(f"\n  Missing odd multiples of 7 (up to {max_H}):")
    for m in range(1, max_H // 7 + 1):
        h = 7 * m
        if h % 2 == 1 and h not in H_counter and h <= max_H:
            note = ""
            if h == 7: note = " = Phi3(2)"
            elif h == 21: note = " = Phi3(4)"
            elif h == 63: note = " = 7*9"
            print(f"    H={h} = 7*{m}{note}")

    # 3. (a1 + 2*a2) mod 7 distribution
    print(f"\n{'='*70}")
    print("(alpha_1 + 2*alpha_2) mod 7 DISTRIBUTION")
    print("=" * 70)
    for s in range(7):
        pct = 100 * sum_mod7[s] / num_samples
        note = " <-- needed for H=0 mod 7" if s == 3 else ""
        print(f"  a1 + 2*a2 mod 7 = {s}: {sum_mod7[s]:7d} ({pct:.2f}%){note}")

    # 4. (a1 mod 7, a2 mod 7) distribution
    print(f"\n{'='*70}")
    print("(alpha_1 mod 7, alpha_2 mod 7) DISTRIBUTION")
    print("=" * 70)
    for a1m in range(7):
        for a2m in range(7):
            cnt = pair_counter.get((a1m, a2m), 0)
            if cnt > 0:
                pct = 100 * cnt / num_samples
                s = (a1m + 2*a2m) % 7
                print(f"  (a1 mod 7, a2 mod 7) = ({a1m},{a2m}): {cnt:6d} ({pct:.1f}%)  sum mod 7 = {s}")

    # 5. All gaps at n=7 (up to 200)
    print(f"\n{'='*70}")
    print("ALL ODD GAPS AT n=7 (up to 200)")
    print("=" * 70)
    for h in range(1, min(201, max_H+1), 2):
        if h not in H_counter:
            note = ""
            if h == 7: note = " = 7*3^0"
            elif h == 21: note = " = 7*3^1"
            elif h == 63: note = " = 7*3^2"
            elif h == 189: note = " = 7*3^3"
            elif h % 7 == 0: note = f" = 7*{h//7}"
            print(f"    H={h}{note}")

    # 6. Distribution of achievable (a1, a2) pairs
    print(f"\n{'='*70}")
    print("MOST COMMON (alpha_1, alpha_2) PAIRS")
    print("=" * 70)
    raw_pairs = Counter()
    for trial in range(50000):  # Quick re-sample for raw pairs
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        a1, a2, _ = compute_alpha_n7(A, n)
        raw_pairs[(a1, a2)] += 1

    for (a1, a2), cnt in sorted(raw_pairs.items(), key=lambda x: -x[1])[:30]:
        H = 1 + 2*a1 + 4*a2
        print(f"  (a1={a1:3d}, a2={a2:3d}) -> H={H:4d}  count={cnt:5d}  H mod 7 = {H%7}")

if __name__ == "__main__":
    main()
