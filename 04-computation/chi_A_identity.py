"""
chi_A = sum (-1)^p |A_p| — alternating sum of allowed path counts.

From euler_char_scaling.py: chi_A = 1 ALWAYS at n=4 (tournament-independent).
At other n, chi_A varies.

Since |A_p| = number of directed Hamiltonian paths on induced (p+1)-subsets,
chi_A = sum_{p=0}^{n-1} (-1)^p * sum_{|S|=p+1} H(T[S])

By switching sums:
chi_A = sum_{k=1}^{n} (-1)^{k-1} * sum_{|S|=k} H(T[S])

For a single vertex set S of size k, H(T[S]) = # Ham paths on T[S].
Define h_k = sum_{|S|=k} H(T[S]).

Then chi_A = sum_{k=1}^{n} (-1)^{k-1} h_k
           = h_1 - h_2 + h_3 - h_4 + ...

For transitive: H(T[S]) = 1 for all S, so h_k = C(n,k).
chi_A(transitive) = sum (-1)^{k-1} C(n,k) = 1 - (1-1)^n = 1 for n >= 1.

For general tournaments: h_k = sum_{|S|=k} H(T[S]).
h_1 = n (always), h_2 = C(n,2) (always, since H(2-tour) = 1).

h_3 = C(n,3) + 2*c3, where c3 = number of directed 3-cycles.
Because H(3-tour) = 1 (transitive) or 3 (cyclic), and 3 = 1+2.

So chi_A = n - C(n,2) + (C(n,3) + 2*c3) - h_4 + ...

The EXCESS over transitive is:
chi_A - 1 = 2*c3 - excess_4 + excess_5 - ...
where excess_k = h_k - C(n,k) = sum_{|S|=k} (H(T[S]) - 1).

Investigate this formula.
"""
import numpy as np
from math import comb
from collections import defaultdict
from itertools import combinations

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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def ham_paths_on_subset(A, subset):
    sub = sorted(subset)
    k = len(sub)
    if k <= 1: return 1
    sub_A = np.zeros((k, k), dtype=int)
    for i, u in enumerate(sub):
        for j, v in enumerate(sub):
            sub_A[i][j] = A[u][v]
    dp = {}
    for v in range(k):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << k):
        bits = bin(mask).count('1')
        if bits == 1: continue
        for v in range(k):
            if not (mask & (1 << v)): continue
            prev_mask = mask ^ (1 << v)
            for u in range(k):
                if (prev_mask & (1 << u)) and sub_A[u][v]:
                    if (mask, v) not in dp:
                        dp[(mask, v)] = 0
                    dp[(mask, v)] += dp.get((prev_mask, u), 0)
    full = (1 << k) - 1
    return sum(dp.get((full, v), 0) for v in range(k))

def main():
    print("=" * 70)
    print("CHI_A IDENTITY: sum (-1)^p |A_p| = sum (-1)^{k-1} h_k")
    print("h_k = sum_{|S|=k} H(T[S])")
    print("=" * 70)

    # Part 1: Compute h_k decomposition for small n
    print("\n--- Part 1: h_k decomposition ---")

    for n in range(3, 8):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        print(f"\n  n={n}:")
        chi_A_vals = []
        excess_k_avgs = defaultdict(float)

        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            h = {}
            excess = {}
            for k in range(1, n+1):
                h_k = sum(ham_paths_on_subset(A, S) for S in combinations(range(n), k))
                h[k] = h_k
                excess[k] = h_k - comb(n, k)

            chi_A = sum((-1)**(k-1) * h[k] for k in range(1, n+1))
            chi_A_vals.append(chi_A)

            for k in range(1, n+1):
                excess_k_avgs[k] += excess[k]

        for k in range(1, n+1):
            avg_excess = excess_k_avgs[k] / N
            c = comb(n, k)
            print(f"    k={k}: C(n,k)={c}, avg excess = {avg_excess:.2f}, "
                  f"avg h_k = {c + avg_excess:.2f}")

        print(f"    chi_A values: min={min(chi_A_vals)}, max={max(chi_A_vals)}, "
              f"all odd: {all(v % 2 == 1 for v in chi_A_vals)}")

        # Check: is chi_A always 1 mod 2?
        print(f"    chi_A mod 2: always {chi_A_vals[0] % 2}: {all(v % 2 == chi_A_vals[0] % 2 for v in chi_A_vals)}")

    # Part 2: chi_A formula via excess
    # chi_A = 1 + sum_{k=1}^{n} (-1)^{k-1} excess_k
    # = 1 + excess_1 - excess_2 + excess_3 - ...
    # excess_1 = 0, excess_2 = 0 always
    # So chi_A = 1 + excess_3 - excess_4 + excess_5 - ...
    # = 1 + 2*c3 - excess_4 + excess_5 - ...
    print("\n" + "=" * 70)
    print("Part 2: chi_A = 1 + 2*c3 - excess_4 + excess_5 - ...")
    print("=" * 70)

    for n in range(3, 8):
        rng = np.random.RandomState(42 + n)
        N = min(200, 2**(n*(n-1)//2))

        print(f"\n  n={n}:")
        for trial in range(min(5, N)):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1

            excesses = {}
            for k in range(1, n+1):
                h_k = sum(ham_paths_on_subset(A, S) for S in combinations(range(n), k))
                excesses[k] = h_k - comb(n, k)

            # Verify excess_3 = 2*c3
            assert excesses[3] == 2 * c3, f"excess_3={excesses[3]} != 2*c3={2*c3}"

            chi_A_computed = 1 + sum((-1)**(k-1) * excesses[k] for k in range(3, n+1))
            chi_A_direct = sum((-1)**(k-1) * (comb(n,k) + excesses[k]) for k in range(1, n+1))
            assert chi_A_computed == chi_A_direct, f"Mismatch: {chi_A_computed} != {chi_A_direct}"

            excess_str = " + ".join(f"({'-' if k % 2 == 0 else '+'}){excesses[k]}" for k in range(3, n+1))
            print(f"    trial {trial}: c3={c3}, excesses={[excesses[k] for k in range(3,n+1)]}, chi_A={chi_A_direct}")

        # Verify excess_3 = 2*c3 for all
        print(f"    Verifying excess_3 = 2*c3 for all {N} samples...", end="")
        ok = True
        for trial in range(N):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)
            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1
            h3 = sum(ham_paths_on_subset(A, S) for S in combinations(range(n), 3))
            if h3 - comb(n, 3) != 2 * c3:
                ok = False
                break
        print(f" {'PASS' if ok else 'FAIL'}")

    # Part 3: Why is chi_A = 1 at n=4?
    # chi_A = 1 + 2*c3 - excess_4
    # excess_4 = H(T) - 1 (since there's only one 4-subset: all vertices)
    # So chi_A = 1 + 2*c3 - (H(T) - 1) = 2 + 2*c3 - H(T)
    # chi_A = 1 means: 2 + 2*c3 - H(T) = 1 => H(T) = 2*c3 + 1
    print("\n" + "=" * 70)
    print("Part 3: At n=4, chi_A = 1 iff H(T) = 2*c3 + 1")
    print("=" * 70)

    n = 4
    total = 2**(n*(n-1)//2)
    all_match = True
    for bits in range(total):
        A = bits_to_adj(bits, n)
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if A[a][b] and A[b][c] and A[c][a]: c3 += 1
            if A[a][c] and A[c][b] and A[b][a]: c3 += 1
        H = ham_paths_on_subset(A, range(n))
        if H != 2 * c3 + 1:
            all_match = False
            print(f"  FAIL at bits={bits}: c3={c3}, H={H}, expected {2*c3+1}")
    if all_match:
        print(f"  VERIFIED: H(T_4) = 2*c3 + 1 for ALL {total} tournaments")
        print(f"  This is equivalent to chi_A = 1 at n=4")
        print(f"  Since c3 in {{0, 1}}: H in {{1, 3}} — but wait, c3=2 gives T_4 cyclic with H=5")

    # Let's check at n=4 what c3 values are possible
    c3_vals = []
    for bits in range(total):
        A = bits_to_adj(bits, n)
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if A[a][b] and A[b][c] and A[c][a]: c3 += 1
            if A[a][c] and A[c][b] and A[b][a]: c3 += 1
        c3_vals.append(c3)
    c3_dist = defaultdict(int)
    for v in c3_vals: c3_dist[v] += 1
    print(f"\n  c3 distribution at n=4: {dict(sorted(c3_dist.items()))}")
    print(f"  H values: {dict(sorted({c3: 2*c3+1 for c3 in c3_dist}.items()))}")

    # Part 4: Does H = 2*c3 + 1 hold at n=5?
    print("\n" + "=" * 70)
    print("Part 4: H vs 2*c3 + 1 at n=5")
    print("=" * 70)

    n = 5
    total = 2**(n*(n-1)//2)
    mismatches = 0
    for bits in range(total):
        A = bits_to_adj(bits, n)
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if A[a][b] and A[b][c] and A[c][a]: c3 += 1
            if A[a][c] and A[c][b] and A[b][a]: c3 += 1
        H = ham_paths_on_subset(A, range(n))
        if H != 2 * c3 + 1:
            mismatches += 1
    print(f"  H(T_5) = 2*c3 + 1: FAILS {mismatches}/{total} times ({100*mismatches/total:.1f}%)")

    # Better: check chi_A formula at n=5
    # chi_A = 1 + 2*c3 - excess_4 + excess_5
    # excess_4 = sum_{|S|=4} (H(T[S]) - 1) = sum of (2*c3(T[S])+1 - 1) = 2 * sum c3(T[S])
    # Wait, this is only true if the n=4 identity holds for each 4-subset...
    # Actually for each 4-subset S, H(T[S]) = 2*c3(S) + 1 by the n=4 identity!
    # So excess_4 = sum_{|S|=4} 2*c3(S) = 2 * C3_4 where C3_4 = sum over 4-subsets of c3

    # Then chi_A = 1 + 2*c3 - 2*C3_4 + (H(T)-1)
    # = H(T) + 2*(c3 - C3_4)

    # Each 3-cycle {a,b,c} is contained in (n-3) 4-subsets.
    # So C3_4 = c3 * (n-3).
    # At n=5: C3_4 = c3 * 2 = 2*c3.
    # chi_A = H(T) + 2*(c3 - 2*c3) = H(T) - 2*c3

    print("\n" + "=" * 70)
    print("Part 5: Recursive formula via n=4 identity")
    print("H(T[S_4]) = 2*c3(S_4) + 1 for ALL 4-subsets S_4")
    print("excess_4 = 2 * sum_{|S|=4} c3(S) = 2 * c3 * (n-3)")
    print("chi_A = 1 + 2*c3 - 2*c3*(n-3) + excess_5 - ...")
    print("=" * 70)

    for n in range(4, 8):
        rng = np.random.RandomState(42 + n)
        N = min(100, 2**(n*(n-1)//2))

        print(f"\n  n={n}:")
        for trial in range(min(10, N)):
            if n <= 5:
                A = bits_to_adj(trial, n)
            else:
                A = random_tournament(n, rng)

            c3 = 0
            for a, b, c in combinations(range(n), 3):
                if A[a][b] and A[b][c] and A[c][a]: c3 += 1
                if A[a][c] and A[c][b] and A[b][a]: c3 += 1

            # Verify: excess_4 = 2 * c3 * (n-3)?
            h4 = sum(ham_paths_on_subset(A, S) for S in combinations(range(n), 4))
            excess_4 = h4 - comb(n, 4)
            predicted_4 = 2 * c3 * (n - 3)

            # Compute chi_A
            chi_A = sum((-1)**(k-1) * sum(ham_paths_on_subset(A, S) for S in combinations(range(n), k)) for k in range(1, n+1))

            print(f"    trial {trial}: c3={c3}, excess_4={excess_4}, predicted={predicted_4}, "
                  f"match={'YES' if excess_4 == predicted_4 else 'NO'}, chi_A={chi_A}")

if __name__ == '__main__':
    main()
