"""
h63_and_moat_structure.py -- kind-pasteur-2026-03-14-S67

The forbidden H values 7 and 21 are Phi_3(KEY_1) and Phi_3(KEY_1^2).
Does the cyclotomic pattern continue? Is H=63 = 7*3^2 also forbidden?
If not, what is the complete "moat structure" beyond H=21?

Exhaustive at n=6 (32768 tournaments), sampling at n=7,8.
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict

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

def compute_H(A, n):
    """Compute H(T) = I(Omega(T), 2) via full odd-cycle decomposition."""
    cycles = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles.append((frozenset(subset), cnt, size))

    # Build independence polynomial coefficients
    # alpha_k = sum of products of k-tuples of cycle counts
    #           over k-tuples of pairwise disjoint cycles
    # H = I(Omega, 2) = sum_k alpha_k * 2^k

    # For efficiency, enumerate independent sets in the conflict graph
    # Using bitmask representation of cycles
    nc = len(cycles)
    if nc == 0:
        return 1

    # Build conflict matrix
    conflict = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i][0] & cycles[j][0]:  # share a vertex
                conflict[i][j] = True
                conflict[j][i] = True

    # Enumerate independent sets via bitmask
    # For each independent set S, contribution = product of (2*cnt_i) for i in S
    # Actually H = 1 + sum over non-empty indep sets of product(2*cnt)
    # More precisely: I(Omega, x) = sum over indep sets product(x * cnt_i)
    # Wait no: I(Omega, x) = sum_{k>=0} alpha_k * x^k
    # where alpha_k = sum over k-indep-sets of product of weights
    # For unweighted: I(G, x) = sum_{k>=0} i_k x^k where i_k = # indep sets of size k
    # For weighted (each vertex has weight = cycle count):
    # I_weighted(G, x) = sum over indep sets S of product_{v in S} (x * w_v)

    # Actually for OCF: H = sum over indep sets S of product_{c in S} (2 * count_c)
    # where count_c = number of directed Hamiltonian cycles on that vertex set.

    H = 1  # empty set contribution
    # Enumerate by BFS/DFS through independent sets
    # For nc up to ~20 this is feasible
    if nc > 25:
        # Use sampling for very large nc
        return -1  # flag as too large

    for mask in range(1, 1 << nc):
        # Check if mask is an independent set
        is_indep = True
        bits = []
        for i in range(nc):
            if mask & (1 << i):
                bits.append(i)
        for ii in range(len(bits)):
            for jj in range(ii+1, len(bits)):
                if conflict[bits[ii]][bits[jj]]:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            prod = 1
            for i in bits:
                prod *= 2 * cycles[i][1]
            H += prod
    return H

def main():
    print("=" * 70)
    print("H=63 ACHIEVABILITY AND MOAT STRUCTURE")
    print("=" * 70)

    # Exhaustive at n=6
    print("\n--- EXHAUSTIVE n=6 ---")
    n = 6
    num_edges = n * (n-1) // 2
    H_counter = Counter()

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
        H = compute_H(A, n)
        H_counter[H] += 1

    max_H_6 = max(H_counter.keys())
    all_H_6 = sorted(H_counter.keys())
    print(f"Total tournaments: {sum(H_counter.values())}")
    print(f"Distinct H values: {len(all_H_6)}")
    print(f"Range: [{min(all_H_6)}, {max_H_6}]")

    # Show all H values up to 70
    print(f"\nAll H values at n=6 (up to 70):")
    for H in range(1, 71):
        if H in H_counter:
            print(f"  H={H:3d}: {H_counter[H]:6d} tournaments", end="")
            if H == 7:
                print("  ← H_forb_1 = Phi_3(2)", end="")
            if H == 21:
                print("  ← H_forb_2 = Phi_3(4)", end="")
            if H == 63:
                print("  ← 7*3^2 (is it forbidden?)", end="")
            print()
        elif H % 2 == 1:  # H is always odd
            if H <= max_H_6:
                print(f"  H={H:3d}: **ABSENT** (gap!)")

    # Find all gaps (missing odd values)
    odd_H_present = set(h for h in H_counter if h % 2 == 1)
    gaps_6 = sorted(h for h in range(1, max_H_6 + 1, 2)
                    if h not in odd_H_present)
    print(f"\nGAPS (missing odd H values at n=6, up to max={max_H_6}):")
    print(f"  {gaps_6}")

    # Now sampling at n=7
    print("\n--- SAMPLING n=7 (all 2097152 tournaments) ---")
    n = 7
    num_edges = n * (n-1) // 2
    H_counter_7 = Counter()
    # n=7 has 2^21 = 2097152 tournaments. Let's do exhaustive.
    # Actually that's a lot with our slow H computation. Sample instead.
    rng = np.random.default_rng(2026_03_14_67)
    num_samples = 500000

    for trial in range(num_samples):
        A = np.zeros((n, n), dtype=np.int8)
        for i in range(n):
            for j in range(i+1, n):
                if rng.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
        H = compute_H(A, n)
        if H >= 0:
            H_counter_7[H] += 1

        if (trial+1) % 100000 == 0:
            print(f"  {trial+1}/{num_samples}...")

    max_H_7 = max(H_counter_7.keys())
    all_H_7 = sorted(H_counter_7.keys())
    print(f"Distinct H values found: {len(all_H_7)}")

    # Check specific values around the forbidden zone
    print(f"\nH values in [1, 30] at n=7:")
    for H in range(1, 31):
        if H in H_counter_7:
            count = H_counter_7[H]
            print(f"  H={H:3d}: {count:6d}", end="")
            if H == 7: print(" ← FORBIDDEN?", end="")
            if H == 21: print(" ← FORBIDDEN?", end="")
            print()
        elif H % 2 == 1:
            print(f"  H={H:3d}: **ABSENT** at n=7 (500k samples)")

    # Check H=63 specifically
    print(f"\nH=63 at n=7: {'FOUND' if 63 in H_counter_7 else 'NOT FOUND'}")
    if 63 in H_counter_7:
        print(f"  Count: {H_counter_7[63]}")

    # Check cyclotomic predictions
    print(f"\n--- CYCLOTOMIC PREDICTIONS ---")
    phi3_vals = [(2, 7), (4, 21), (8, 73)]
    for x, phi3_x in phi3_vals:
        in_6 = phi3_x in H_counter
        in_7 = phi3_x in H_counter_7
        print(f"  Phi_3({x}) = {phi3_x}: n=6 {'PRESENT' if in_6 else 'ABSENT'}, "
              f"n=7 {'PRESENT' if in_7 else 'ABSENT (500k)'}")

    # Check 7*3^k sequence
    print(f"\n--- SEQUENCE H = 7·3^k ---")
    for k in range(6):
        H = 7 * 3**k
        in_6 = H in H_counter
        in_7 = H in H_counter_7
        status_6 = f"{H_counter.get(H, 0):6d}" if H <= max_H_6 else "too big"
        status_7 = f"{H_counter_7.get(H, 0):6d}" if H <= max_H_7 else "too big"
        print(f"  k={k}: H = 7·3^{k} = {H:6d}  |  n=6: {status_6}  |  n=7: {status_7}")

    # Additional: what are the "permanent gaps" (absent at both n=6 and n=7)?
    permanent_gaps = []
    for H in range(1, min(max_H_6, 100) + 1, 2):
        if H not in H_counter and H not in H_counter_7:
            permanent_gaps.append(H)
    print(f"\nPERMANENT GAPS (absent at BOTH n=6 and n=7, H<100):")
    print(f"  {permanent_gaps}")

    # Deeper analysis: check ALL values up to 100 at n=7
    print(f"\nFull gap analysis at n=7 (H values in [1,100]):")
    for H in range(1, 101, 2):
        if H not in H_counter_7 and H <= max_H_7:
            note = ""
            if H == 7: note = " = Phi_3(2)"
            elif H == 21: note = " = Phi_3(4)"
            elif H == 63: note = " = 7·3²"
            elif H == 73: note = " = Phi_3(8)"
            print(f"  H={H:3d}: ABSENT{note}")

    print(f"\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)

if __name__ == "__main__":
    main()
