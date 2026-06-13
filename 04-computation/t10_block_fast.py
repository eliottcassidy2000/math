"""
t10_block_fast.py -- kind-pasteur-2026-03-14-S66

Fast check: does alpha_1=8 ever have alpha_2=1?
Does alpha_1=6 ever have alpha_2=2?
These are the 2 remaining empirical cases for the T=10 universal block.
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

def main():
    print("=" * 70)
    print("FAST T=10 BLOCK: alpha_1=8 and alpha_1=6 forcing at n=7,8")
    print("=" * 70)

    for n in [7, 8]:
        print(f"\n{'='*50}")
        print(f"n = {n}")
        print(f"{'='*50}")

        rng = np.random.default_rng(2026_03_14_50 + n)

        a8_a2 = Counter()
        a6_a2 = Counter()
        t10_found = 0

        num_samples = {7: 300000, 8: 200000}[n]

        for trial in range(num_samples):
            A = np.zeros((n, n), dtype=np.int8)
            for i in range(n):
                for j in range(i+1, n):
                    if rng.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1

            # Quick c3
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or \
                   (A[a][c] and A[c][b] and A[b][a]):
                    c3 += 1

            # alpha_1=8 needs c3<=8 (could have 5-cycle contributions)
            # alpha_1=6 needs c3<=6
            if c3 > 10:
                continue

            cycles = []
            for size in range(3, n+1, 2):
                for subset in combinations(range(n), size):
                    cnt = count_directed_hamcycles(A, list(subset))
                    if cnt > 0:
                        cycles.append((frozenset(subset), cnt, size))

            alpha_1 = sum(cnt for _, cnt, _ in cycles)

            if alpha_1 in (6, 8):
                alpha_2 = 0
                for i in range(len(cycles)):
                    for j in range(i+1, len(cycles)):
                        if len(cycles[i][0] & cycles[j][0]) == 0:
                            alpha_2 += cycles[i][1] * cycles[j][1]

                if alpha_1 == 8:
                    a8_a2[alpha_2] += 1
                else:
                    a6_a2[alpha_2] += 1

                T = alpha_1 + 2 * alpha_2
                if T == 10:
                    t10_found += 1
                    scores = tuple(sorted(int(sum(A[i])) for i in range(n)))
                    print(f"  *** T=10 FOUND! a1={alpha_1}, a2={alpha_2}, scores={scores}")

            if (trial+1) % 100000 == 0:
                print(f"  {trial+1}/{num_samples}...")

        print(f"\nT=10 found at n={n}: {t10_found}")

        total_a8 = sum(a8_a2.values())
        print(f"\nalpha_1=8 at n={n}: {total_a8} found")
        if total_a8:
            for a2 in sorted(a8_a2.keys()):
                T = 8 + 2*a2
                print(f"  alpha_2={a2}: {a8_a2[a2]} ({100*a8_a2[a2]/total_a8:.1f}%), T={T}")

        total_a6 = sum(a6_a2.values())
        print(f"\nalpha_1=6 at n={n}: {total_a6} found")
        if total_a6:
            for a2 in sorted(a6_a2.keys()):
                T = 6 + 2*a2
                print(f"  alpha_2={a2}: {a6_a2[a2]} ({100*a6_a2[a2]/total_a6:.1f}%), T={T}")

    # Also check n=9 with fewer samples
    print(f"\n{'='*50}")
    print(f"n = 9 (light sampling)")
    print(f"{'='*50}")
    n = 9
    rng = np.random.default_rng(2026_03_14_59)
    a8_a2 = Counter()
    a6_a2 = Counter()
    t10_found = 0

    for trial in range(100000):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or \
               (A[a][c] and A[c][b] and A[b][a]):
                c3 += 1
        if c3 > 10:
            continue

        cycles = []
        for size in range(3, min(n+1, 10), 2):  # cap at 9
            for subset in combinations(range(n), size):
                cnt = count_directed_hamcycles(A, list(subset))
                if cnt > 0:
                    cycles.append((frozenset(subset), cnt, size))

        alpha_1 = sum(cnt for _, cnt, _ in cycles)

        if alpha_1 in (6, 8):
            alpha_2 = 0
            for i in range(len(cycles)):
                for j in range(i+1, len(cycles)):
                    if len(cycles[i][0] & cycles[j][0]) == 0:
                        alpha_2 += cycles[i][1] * cycles[j][1]

            if alpha_1 == 8:
                a8_a2[alpha_2] += 1
            else:
                a6_a2[alpha_2] += 1

            T = alpha_1 + 2 * alpha_2
            if T == 10:
                t10_found += 1
                print(f"  *** T=10 at n=9! a1={alpha_1}, a2={alpha_2}")

        if (trial+1) % 50000 == 0:
            print(f"  {trial+1}/100000...")

    print(f"\nT=10 found at n=9: {t10_found}")
    total_a8 = sum(a8_a2.values())
    print(f"alpha_1=8 at n=9: {total_a8}")
    if total_a8:
        for a2 in sorted(a8_a2.keys()):
            print(f"  alpha_2={a2}: {a8_a2[a2]}, T={8+2*a2}")
    total_a6 = sum(a6_a2.values())
    print(f"alpha_1=6 at n=9: {total_a6}")
    if total_a6:
        for a2 in sorted(a6_a2.keys()):
            print(f"  alpha_2={a2}: {a6_a2[a2]}, T={6+2*a2}")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)
    print()
    print("For H=21 (T=10) to be achievable, one of these must hold:")
    print("  (10,0): PROVED IMPOSSIBLE (alpha_1=10 forces alpha_2>=2)")
    print("  (8,1):  See alpha_1=8 data above")
    print("  (6,2):  See alpha_1=6 data above")
    print("  (4,3):  PROVED IMPOSSIBLE (alpha_1=4 has alpha_2 in {0,4} only)")
    print("  (2,4):  PROVED IMPOSSIBLE (2 cycles => alpha_2<=1)")
    print("  (0,5):  PROVED IMPOSSIBLE (0 cycles => alpha_2=0)")

if __name__ == "__main__":
    main()
