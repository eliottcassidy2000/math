"""
t10_universal_block.py -- kind-pasteur-2026-03-14-S66

T=10 Universal Block Theorem (HYP-1048):
  T = alpha_1 + 2*alpha_2 = 10 is NEVER achieved by any tournament at any n.

Exhaustive verification at n=6 shows ALL 6 decompositions are blocked:
  (10,0): alpha_1=10 forces alpha_2=2 (PROVED)
  (8,1):  alpha_1=8 only has alpha_2=0 at n=6
  (6,2):  alpha_1=6 only has alpha_2 in {0,1} at n=6
  (4,3):  alpha_1=4 only has alpha_2=0 at n=6 (alpha_2 in {0,4} at n=7-9)
  (2,4):  PROVED impossible (2 cycles can have at most 1 disjoint pair)
  (0,5):  PROVED impossible (0 cycles => alpha_2=0)

Check: do (8,1) and (6,2) remain blocked at n=7,8?
Also: what alpha_2 values does alpha_1=8 and alpha_1=6 achieve at n=7,8?
"""

import numpy as np
from itertools import combinations
from collections import Counter

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

def analyze_tournament(A, n):
    """Return (alpha_1, alpha_2, alpha_3, cycle_list)."""
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

    alpha_3 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) > 0:
                continue
            for k_idx in range(j+1, len(cycles)):
                if len(cycles[i][0] & cycles[k_idx][0]) == 0 and \
                   len(cycles[j][0] & cycles[k_idx][0]) == 0:
                    alpha_3 += cycles[i][1] * cycles[j][1] * cycles[k_idx][1]

    return alpha_1, alpha_2, alpha_3, cycles

def main():
    print("=" * 70)
    print("T=10 UNIVERSAL BLOCK: CHECKING (8,1) AND (6,2) AT n=7,8")
    print("=" * 70)

    for n in [7, 8]:
        print(f"\n{'='*50}")
        print(f"n = {n}")
        print(f"{'='*50}")

        rng = np.random.default_rng(2026_03_14_10 + n)

        # Track alpha_2 distributions for specific alpha_1 values
        a1_to_a2 = {a1: Counter() for a1 in [4, 6, 8, 10]}
        t10_found = 0

        num_samples = {7: 1000000, 8: 500000}[n]

        for trial in range(num_samples):
            A = np.zeros((n, n), dtype=np.int8)
            for i in range(n):
                for j in range(i+1, n):
                    if rng.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1

            # Quick c3 count for filtering
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    c3 += 1

            # For alpha_1 in {4,6,8,10}, c3 gives a rough bound
            # c3 <= 12 covers all these cases
            if c3 > 12:
                continue

            alpha_1, alpha_2, alpha_3, cycles = analyze_tournament(A, n)

            if alpha_1 in a1_to_a2:
                a1_to_a2[alpha_1][alpha_2] += 1

            T = alpha_1 + 2 * alpha_2
            if T == 10:
                t10_found += 1
                print(f"  *** T=10 FOUND! alpha_1={alpha_1}, alpha_2={alpha_2}, alpha_3={alpha_3}")
                scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
                print(f"      scores={scores}, c3={c3}")
                for vs, cnt, sz in cycles:
                    print(f"      Cycle: {sorted(vs)} (size {sz}, {cnt} dir)")

            if (trial+1) % 200000 == 0:
                print(f"  {trial+1}/{num_samples}...")

        print(f"\nResults at n={n} ({num_samples} samples):")
        print(f"T=10 found: {t10_found}")

        for a1 in sorted(a1_to_a2.keys()):
            a2_dist = a1_to_a2[a1]
            total = sum(a2_dist.values())
            if total > 0:
                print(f"\nalpha_1={a1}: {total} found")
                for a2 in sorted(a2_dist.keys()):
                    pct = 100 * a2_dist[a2] / total
                    T = a1 + 2*a2
                    H_min = 1 + 2*a1 + 4*a2
                    print(f"  alpha_2={a2}: {a2_dist[a2]} ({pct:.1f}%), T={T}, H_min={H_min}")
            else:
                print(f"\nalpha_1={a1}: 0 found")

    # Theoretical summary
    print(f"\n{'='*70}")
    print("SUMMARY: T=10 BLOCK STATUS")
    print("=" * 70)
    print()
    print("Decomposition  | n=6 (exh.) | n=7 (1M)  | n=8 (500k) | Status")
    print("-" * 70)
    print("(10,0)         | alpha_2=2  | alpha_2=2 | alpha_2=2  | PROVED")
    print("(8,1)          | alpha_2=0  | see above | see above  | ???")
    print("(6,2)          | alpha_2={0,1}| see above| see above  | ???")
    print("(4,3)          | alpha_2=0  | alpha_2=0 | alpha_2={0,4}| STRONG")
    print("(2,4)          | IMPOSSIBLE | IMPOSSIBLE| IMPOSSIBLE | PROVED")
    print("(0,5)          | IMPOSSIBLE | IMPOSSIBLE| IMPOSSIBLE | PROVED")

if __name__ == "__main__":
    main()
