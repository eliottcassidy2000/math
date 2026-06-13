"""
h63_fast.py -- kind-pasteur-2026-03-14-S67

Fast check: is H=63 achievable? What about other "predicted forbidden" values?
Uses fast cycle counting instead of exhaustive independence polynomial.

At n=6 exhaustively, at n=7,8 by sampling.
H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
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

def compute_H_fast(A, n):
    """Compute H via alpha_k for k=1,2,3."""
    cycles = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt, size))

    alpha_1 = sum(cnt for _, cnt, _ in cycles)

    # alpha_2 = weighted independent pairs
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) == 0:
                alpha_2 += cycles[i][1] * cycles[j][1]

    # alpha_3 = weighted independent triples
    alpha_3 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i][0] & cycles[j][0]) > 0:
                continue
            for k in range(j+1, len(cycles)):
                if len(cycles[i][0] & cycles[k][0]) == 0 and \
                   len(cycles[j][0] & cycles[k][0]) == 0:
                    alpha_3 += cycles[i][1] * cycles[j][1] * cycles[k][1]

    H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
    return H, alpha_1, alpha_2, alpha_3

def main():
    print("=" * 70)
    print("H=63 AND MOAT STRUCTURE — FAST VERSION")
    print("=" * 70)

    # n=6 exhaustive (only 3-cycles and 5-cycles exist)
    print("\n--- EXHAUSTIVE n=6 ---")
    n = 6
    num_edges = n * (n - 1) // 2  # 15
    H_counter = Counter()
    alpha_data = {}

    for bits in range(1 << num_edges):
        A = np.zeros((n, n), dtype=np.int8)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        H, a1, a2, a3 = compute_H_fast(A, n)
        H_counter[H] += 1

        if bits % 8192 == 0:
            pct = 100 * bits / (1 << num_edges)
            print(f"  {pct:.0f}%...", flush=True)

    max_H = max(H_counter.keys())
    print(f"\nTotal tournaments: {sum(H_counter.values())}")
    print(f"Distinct H values: {len(H_counter)}")
    print(f"Min H: {min(H_counter.keys())}, Max H: {max_H}")

    # All H values
    print(f"\nAll H values at n=6:")
    for H in sorted(H_counter.keys()):
        note = ""
        if H == 7: note = " <-- Phi3(2) FORBIDDEN?"
        elif H == 21: note = " <-- Phi3(4) FORBIDDEN?"
        elif H == 63: note = " <-- 7*3^2 FORBIDDEN?"
        elif H == 73: note = " <-- Phi3(8)"
        print(f"  H={H:4d}: {H_counter[H]:6d}{note}")

    # Gap analysis
    all_H = sorted(H_counter.keys())
    print(f"\nGAPS (odd values NOT achieved, up to {max_H}):")
    gaps = []
    for h in range(1, max_H + 1, 2):
        if h not in H_counter:
            gaps.append(h)
    print(f"  {gaps}")

    # n=7 sampling
    print(f"\n--- SAMPLING n=7 (200k) ---")
    n = 7
    rng = np.random.default_rng(2026_03_14_97)
    H_counter_7 = Counter()

    for trial in range(200000):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        H, a1, a2, a3 = compute_H_fast(A, n)
        H_counter_7[H] += 1

        if (trial+1) % 50000 == 0:
            print(f"  {trial+1}/200000...")

    max_H_7 = max(H_counter_7.keys())
    print(f"Distinct H values: {len(H_counter_7)}")

    # Check specific values
    print(f"\nKey H values at n=7:")
    for target in [7, 21, 63, 73]:
        found = target in H_counter_7
        count = H_counter_7.get(target, 0)
        print(f"  H={target}: {'FOUND' if found else 'NOT FOUND'} ({count} times)")

    print(f"\nAll ABSENT odd H values <= 100 at n=7:")
    for h in range(1, 101, 2):
        if h not in H_counter_7 and h <= max_H_7:
            note = ""
            if h == 7: note = " = Phi3(2)"
            elif h == 21: note = " = Phi3(4)"
            elif h == 63: note = " = 7*3^2"
            elif h == 73: note = " = Phi3(8)"
            print(f"  H={h}{note}")

    # Summary
    print(f"\n{'='*70}")
    print("MOAT STRUCTURE SUMMARY")
    print("=" * 70)
    print(f"\nH values present at n=6: {sorted(H_counter.keys())}")
    print(f"\nPERMANENT GAPS (absent at both n=6 and n=7, H <= 100):")
    for h in range(1, 101, 2):
        if h not in H_counter and h not in H_counter_7:
            if h <= max_H_7:
                note = ""
                if h == 7: note = " = Phi3(KEY1) = I(K3,2)"
                elif h == 21: note = " = Phi3(KEY1^2) = I(K3+K1,2)"
                print(f"  H={h}{note}")

if __name__ == "__main__":
    main()
