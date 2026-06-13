"""
moat_n8_analysis.py -- kind-pasteur-2026-03-14-S67

Targeted n=8 moat analysis:
1. Does H=63 remain forbidden at n=8? (Expected: YES)
2. Is H=189 achievable at n=8? (Expected: YES, since H(T_7)=189)
3. Is H=567 = 7*3^4 forbidden at n=8?
4. Full gap structure of achievable H values at n=8

Uses fast alpha computation: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + 16*alpha_4
At n=8: max cycle size = 7 (odd), so alpha_4 = 0 for most tournaments.
Actually at n=8: odd cycles on subsets of size 3, 5, 7.
Four disjoint 3-cycles impossible (need 12 > 8 vertices).
Max independent set: two disjoint 3-cycles (6 vertices) + nothing else.
Or one 3-cycle + one 5-cycle on disjoint vertices = 8 vertices exactly.
So alpha_3 is possible, alpha_4 is impossible at n=8 (would need >= 12 vertices).
Wait: alpha_k counts independent sets of size k in Omega(T).
A 7-cycle uses all but 1 vertex. Can't be disjoint from anything.
Max independent set size in Omega(T) at n=8:
- Two disjoint 3-cycles = 6 vertices => size 2.
- One 3-cycle + one 5-cycle on disjoint vertices = 8 vertices => size 2.
- Three disjoint 3-cycles would need 9 > 8 vertices => impossible.
So alpha_3 = 0 at n=8! H = 1 + 2*alpha_1 + 4*alpha_2.

Wait no. alpha_k counts independent sets of size k in Omega(T), where
Omega(T) is the conflict graph. Each vertex of Omega is a directed odd cycle.
An independent set is a collection of pairwise vertex-disjoint cycles.
At n=8, we can have:
- Single cycle: 3,5,7-cycle (alpha_1 counts total weighted)
- Pair of disjoint cycles: (3,3), (3,5) -- both fit in 8 vertices
- Triple of disjoint cycles: (3,3,...) would need a third on remaining 2 vertices = impossible
So max IS size = 2. alpha_3 = 0 at n=8.

But wait: alpha_1 counts the total number of directed cycles (weighted),
and alpha_2 counts pairs of vertex-disjoint directed cycles.
At n=8, the max is pairs: (3,3) using 6 of 8, or (3,5) using all 8.

So H = 1 + 2*alpha_1 + 4*alpha_2 at n=8.
"""

import numpy as np
from itertools import combinations
from collections import Counter
import time

def count_directed_hamcycles(A, vertices):
    """Count directed Hamiltonian cycles on a subset of vertices."""
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

def compute_H_n8(A, n):
    """Compute H = 1 + 2*alpha_1 + 4*alpha_2 for n=8."""
    # Enumerate all directed odd cycles on subsets
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

    H = 1 + 2*alpha_1 + 4*alpha_2
    return H, alpha_1, alpha_2

def main():
    print("=" * 70)
    print("N=8 MOAT ANALYSIS — FORBIDDEN 7*3^k SEQUENCE")
    print("=" * 70)

    n = 8
    rng = np.random.default_rng(2026_03_14_67)

    # Track key values
    H_counter = Counter()
    targets = {7, 21, 63, 189, 567}
    target_hits = {t: 0 for t in targets}

    num_samples = 100000
    start = time.time()

    for trial in range(num_samples):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        H, a1, a2 = compute_H_n8(A, n)
        H_counter[H] += 1

        for t in targets:
            if H == t:
                target_hits[t] += 1

        if (trial+1) % 10000 == 0:
            elapsed = time.time() - start
            rate = (trial+1) / elapsed
            print(f"  {trial+1}/{num_samples} ({rate:.0f}/sec)...", flush=True)

    elapsed = time.time() - start
    print(f"\nCompleted {num_samples} samples in {elapsed:.1f}s")

    max_H = max(H_counter.keys())
    min_H = min(H_counter.keys())
    print(f"Distinct H values: {len(H_counter)}")
    print(f"Min H: {min_H}, Max H: {max_H}")

    # Key target check
    print(f"\n{'='*70}")
    print("FORBIDDEN SEQUENCE 7*3^k CHECK")
    print("=" * 70)
    for k in range(5):
        val = 7 * 3**k
        found = target_hits.get(val, H_counter.get(val, 0))
        status = "FOUND" if found > 0 else "NOT FOUND"
        print(f"  7*3^{k} = {val:5d}: {status} ({found} hits)")

    # Full gap analysis
    print(f"\n{'='*70}")
    print("GAP ANALYSIS (odd H values NOT achieved)")
    print("=" * 70)
    all_gaps = []
    for h in range(1, min(max_H + 1, 601), 2):
        if h not in H_counter:
            all_gaps.append(h)

    print(f"Total gaps (odd, <= {min(max_H, 600)}): {len(all_gaps)}")
    print(f"\nGaps <= 200:")
    for h in all_gaps:
        if h <= 200:
            note = ""
            if h == 7: note = " = 7*3^0 = Phi3(2)"
            elif h == 21: note = " = 7*3^1 = Phi3(4)"
            elif h == 63: note = " = 7*3^2"
            elif h == 189: note = " = 7*3^3 = H(T_7)"
            elif h == 567: note = " = 7*3^4"
            elif h % 7 == 0: note = f" = 7*{h//7}"
            print(f"    H={h}{note}")

    print(f"\nGaps 200-600:")
    for h in all_gaps:
        if 200 < h <= 600:
            note = ""
            if h == 567: note = " = 7*3^4"
            elif h % 7 == 0: note = f" = 7*{h//7}"
            print(f"    H={h}{note}")

    # Check 7*3^k pattern specifically
    print(f"\n{'='*70}")
    print("7*3^k MOAT BOUNDARY")
    print("=" * 70)
    for k in range(6):
        val = 7 * 3**k
        if val > max_H:
            print(f"  7*3^{k} = {val}: beyond max H={max_H} (cannot test)")
            break
        found = val in H_counter
        print(f"  7*3^{k} = {val}: {'ACHIEVABLE' if found else 'FORBIDDEN'}")

    # Distribution around gaps
    print(f"\n{'='*70}")
    print("DISTRIBUTION NEAR FORBIDDEN VALUES")
    print("=" * 70)
    for target in [7, 21, 63, 189]:
        if target <= max_H + 10:
            print(f"\n  Near H={target}:")
            for h in range(max(1, target-6), min(max_H+1, target+8), 2):
                count = H_counter.get(h, 0)
                bar = "#" * min(count // 100, 40)
                mark = " <<<" if h == target else ""
                print(f"    H={h:4d}: {count:6d} {bar}{mark}")

    # H mod 7 distribution
    print(f"\n{'='*70}")
    print("H mod 7 DISTRIBUTION")
    print("=" * 70)
    mod7 = Counter()
    for h, cnt in H_counter.items():
        mod7[h % 7] += cnt
    for r in range(7):
        pct = 100 * mod7[r] / num_samples
        print(f"  H ≡ {r} mod 7: {mod7[r]:6d} ({pct:.1f}%)")

    # H mod 21 distribution
    print(f"\n  H mod 21 for multiples of 7:")
    mod21 = Counter()
    for h, cnt in H_counter.items():
        if h % 7 == 0:
            mod21[h % 21] += cnt
    for r in sorted(mod21.keys()):
        print(f"    H ≡ {r} mod 21: {mod21[r]:6d}")

if __name__ == "__main__":
    main()
