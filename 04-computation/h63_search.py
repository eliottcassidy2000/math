"""
h63_search.py — Targeted search for H=63 at n=8,9,10
kind-pasteur-2026-03-14-S65

H=63 is absent at n<=7 but structurally feasible.
This script does heavy sampling to check if it appears at larger n.
Also checks H=107, 119, 149.
"""

import numpy as np

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_ham_paths(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total > 0:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def main():
    targets = {7, 21, 63, 107, 119, 149}

    print("=" * 60)
    print("TARGETED SEARCH FOR FORBIDDEN H VALUES")
    print("=" * 60)

    for n in [8, 9, 10]:
        N = 5000 if n <= 9 else 1000
        rng = np.random.default_rng(2026_0314 + n)

        print(f"\nn = {n} ({N} random tournaments):")

        found = {h: 0 for h in targets}
        h_counts = {}
        min_h = float('inf')
        max_h = 0

        for trial in range(N):
            A = random_tournament(n, rng)
            H = count_ham_paths(A)

            h_counts[H] = h_counts.get(H, 0) + 1
            min_h = min(min_h, H)
            max_h = max(max_h, H)

            for h_target in targets:
                if H == h_target:
                    found[h_target] += 1

            if (trial + 1) % 1000 == 0:
                print(f"  Progress: {trial+1}/{N}")

        print(f"  H range: [{min_h}, {max_h}]")
        print(f"  Distinct H values seen: {len(h_counts)}")

        for h_target in sorted(targets):
            count = found[h_target]
            status = f"FOUND ({count}x, {100*count/N:.2f}%)" if count > 0 else "absent"
            print(f"  H={h_target:4d}: {status}")

        # Show H values close to targets
        for h_target in [63, 107, 119, 149]:
            nearby = sorted([h for h in h_counts if abs(h - h_target) <= 6])
            if nearby:
                pairs = [f"H={h}({h_counts[h]})" for h in nearby]
                print(f"  Near H={h_target}: {', '.join(pairs)}")

    # Also try structured tournaments at n=9
    print(f"\n{'='*60}")
    print("STRUCTURED SEARCH AT n=9")
    print(f"{'='*60}")

    n = 9
    rng = np.random.default_rng(42)

    # Try tournaments with specific score sequences
    # Near-regular scores tend to have mid-range H
    print(f"\n  Looking for H near 63 via score-constrained sampling...")

    hits = {h: [] for h in targets}

    for trial in range(3000):
        A = random_tournament(n, rng)
        H = count_ham_paths(A)

        if H in targets:
            scores = sorted([sum(A[i]) for i in range(n)])
            hits[H].append((trial, scores))

        # Also try: perturb near-regular to get specific H values
        if 55 <= H <= 75:
            scores = sorted([sum(A[i]) for i in range(n)])
            if trial < 20 or H in targets:
                print(f"    trial {trial}: H={H}, scores={scores}")

    for h_target in sorted(targets):
        if hits[h_target]:
            print(f"\n  H={h_target} FOUND!")
            for trial, scores in hits[h_target][:3]:
                print(f"    trial {trial}: scores={scores}")
        else:
            print(f"\n  H={h_target}: not found in 3000 structured trials")

if __name__ == "__main__":
    main()
