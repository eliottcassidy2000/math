"""
vassiliev_tournament.py -- kind-pasteur-2026-03-14-S69
Tournament Vassiliev-type (finite-type) invariants.

In knot theory: a Vassiliev invariant of type k is one that vanishes
on all k-singular knots (knots with k double-point self-intersections).

Tournament analogy: a "singular arc" is one where we average over both
orientations. An invariant v is of finite type k if for every tournament T
and every set S of k+1 arcs: sum_{S' subset S} (-1)^|S'| v(T_{flip S'}) = 0
where T_{flip S'} flips the arcs in S'.

This is exactly the discrete derivative (finite difference):
Delta_S v(T) = sum over subsets S' of S: (-1)^|S'| v(T xor S') = 0

KEY: H(T) is invariant of type... what? If type 1, then every single arc flip
gives Delta_e H = H(T+) - H(T-) that depends only on the "resolved" version.
But Delta_e H varies, so H is NOT type 0.

Is there a k such that Delta_{S} H = 0 for all |S| = k+1?

Also explore: the "spectral" Vassiliev invariants from the transfer matrix.
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

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

def compute_H_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def arcs_list(n):
    """List all C(n,2) arcs as (i,j) with i<j."""
    return [(i, j) for i in range(n) for j in range(i+1, n)]

def flip_arcs(A, arc_set):
    """Flip the arcs in arc_set."""
    B = A.copy()
    for (i, j) in arc_set:
        B[i][j], B[j][i] = B[j][i], B[i][j]
    return B

def finite_diff(A, n, arc_set, invariant_fn):
    """
    Compute the k-th order finite difference of invariant_fn at T = A,
    with respect to the arc set arc_set.
    Delta_S v(T) = sum_{S' subset S} (-1)^|S'| v(T_{flip S'})
    """
    k = len(arc_set)
    total = 0
    for mask in range(1 << k):
        subset = [arc_set[i] for i in range(k) if mask & (1 << i)]
        B = flip_arcs(A, subset)
        sign = (-1) ** bin(mask).count('1')
        total += sign * invariant_fn(B, n)
    return total

def main():
    print("=" * 70)
    print("VASSILIEV-TYPE INVARIANTS FOR TOURNAMENTS")
    print("kind-pasteur-2026-03-14-S69")
    print("=" * 70)

    # ====================================
    # PART 1: Is H(T) of finite type?
    # ====================================
    print("\n" + "=" * 70)
    print("PART 1: FINITE-TYPE ORDER OF H(T)")
    print("  Check: for each k, does Delta_S H = 0 for ALL |S|=k+1?")
    print("=" * 70)

    for n in [3, 4, 5]:
        print(f"\n--- n = {n} ---")
        arcs = arcs_list(n)
        m = len(arcs)

        # For each order k=1,2,...,m, check if Delta_S H = 0 for all S of size k
        for k in range(1, min(m + 1, 7)):
            nonzero_count = 0
            total_tests = 0
            max_abs_delta = 0

            # Sample arc subsets (or exhaustive for small)
            arc_subsets = list(combinations(range(m), k))
            if len(arc_subsets) > 200:
                # Sample
                import random
                random.seed(42)
                arc_subsets = random.sample(arc_subsets, 200)

            for bits in range(min(2**m, 100)):
                A = bits_to_adj(bits, n)
                for subset_idx in arc_subsets[:50]:
                    arc_subset = [arcs[i] for i in subset_idx]
                    delta = finite_diff(A, n, arc_subset, compute_H_dp)
                    total_tests += 1
                    if delta != 0:
                        nonzero_count += 1
                    max_abs_delta = max(max_abs_delta, abs(delta))

            if nonzero_count == 0:
                print(f"  k={k}: ALL ZERO ({total_tests} tests). H is type {k-1}!")
                break
            else:
                print(f"  k={k}: {nonzero_count}/{total_tests} nonzero, max|Delta|={max_abs_delta}")

    # ====================================
    # PART 2: The k=1 finite difference = arc-flip delta
    # ====================================
    print("\n" + "=" * 70)
    print("PART 2: k=1 FINITE DIFFERENCE (ARC-FLIP DELTAS)")
    print("  Delta_e H = H(T) - H(T flip e)")
    print("  This is our familiar arc-flip spectrum")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")
        arcs = arcs_list(n)
        delta_distribution = Counter()

        for bits in range(2**(n*(n-1)//2)):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            for arc in arcs:
                B = flip_arcs(A, [arc])
                H_flip = compute_H_dp(B, n)
                delta = H - H_flip
                delta_distribution[delta] += 1

        print(f"  Delta_e H distribution: {dict(sorted(delta_distribution.items()))}")
        print(f"  All deltas even: {all(d % 2 == 0 for d in delta_distribution.keys())}")
        print(f"  Symmetric: {all(delta_distribution.get(d,0) == delta_distribution.get(-d,0) for d in delta_distribution)}")

    # ====================================
    # PART 3: The k=2 finite difference
    # ====================================
    print("\n" + "=" * 70)
    print("PART 3: k=2 FINITE DIFFERENCE")
    print("  Delta_{e1,e2} H = H(T) - H(T^e1) - H(T^e2) + H(T^{e1,e2})")
    print("  If all zero, H would be of type 1 (linear in arc flips)")
    print("=" * 70)

    for n in [4, 5]:
        print(f"\n--- n = {n} ---")
        arcs = arcs_list(n)
        m = len(arcs)
        delta2_distribution = Counter()
        total_tests = 0

        for bits in range(min(2**(n*(n-1)//2), 200)):
            A = bits_to_adj(bits, n)
            for i in range(m):
                for j in range(i+1, m):
                    arc_pair = [arcs[i], arcs[j]]
                    delta2 = finite_diff(A, n, arc_pair, compute_H_dp)
                    delta2_distribution[delta2] += 1
                    total_tests += 1

        print(f"  Total tests: {total_tests}")
        print(f"  Delta_2 distribution: {dict(sorted(delta2_distribution.items()))}")
        print(f"  All zero: {all(d == 0 for d in delta2_distribution.keys())}")

        # Check: are adjacent vs non-adjacent arc pairs different?
        print(f"\n  Adjacent vs non-adjacent arc pairs:")
        adj_delta = Counter()
        nonadj_delta = Counter()

        for bits in range(min(2**(n*(n-1)//2), 100)):
            A = bits_to_adj(bits, n)
            for i in range(m):
                for j in range(i+1, m):
                    a1 = arcs[i]
                    a2 = arcs[j]
                    shared = len(set(a1) & set(a2))
                    arc_pair = [a1, a2]
                    delta2 = finite_diff(A, n, arc_pair, compute_H_dp)

                    if shared > 0:
                        adj_delta[delta2] += 1
                    else:
                        nonadj_delta[delta2] += 1

        print(f"    Adjacent arcs (shared vertex): {dict(sorted(adj_delta.items()))}")
        print(f"    Non-adjacent arcs (disjoint): {dict(sorted(nonadj_delta.items()))}")

    # ====================================
    # PART 4: Disjoint arc k=2 differences
    # ====================================
    print("\n" + "=" * 70)
    print("PART 4: DISJOINT ARC PAIRS — IS H ADDITIVE ON DISJOINT FLIPS?")
    print("  If Delta_{e1,e2} H = 0 for disjoint arcs, then H decomposes")
    print("  as a sum of 'local' contributions from vertex neighborhoods")
    print("=" * 70)

    for n in [5, 6]:
        print(f"\n--- n = {n} ---")
        arcs = arcs_list(n)
        m = len(arcs)

        disjoint_delta2 = Counter()
        total = 0

        for bits_idx in range(min(2**(n*(n-1)//2), 500)):
            if n >= 6:
                import random
                random.seed(bits_idx)
                bits = random.randint(0, 2**(n*(n-1)//2) - 1)
            else:
                bits = bits_idx
            A = bits_to_adj(bits, n)

            # Find disjoint arc pairs
            for i in range(m):
                for j in range(i+1, m):
                    a1 = arcs[i]
                    a2 = arcs[j]
                    if len(set(a1) & set(a2)) == 0:  # disjoint
                        delta2 = finite_diff(A, n, [a1, a2], compute_H_dp)
                        disjoint_delta2[delta2] += 1
                        total += 1

            if total > 10000:
                break

        print(f"  Total disjoint pair tests: {total}")
        print(f"  Disjoint Delta_2 distribution: {dict(sorted(disjoint_delta2.items()))}")
        all_zero = all(d == 0 for d in disjoint_delta2.keys())
        print(f"  All zero (additive on disjoint flips): {all_zero}")

    # ====================================
    # PART 5: 3-cycle reversal = 3-arc simultaneous flip
    # ====================================
    print("\n" + "=" * 70)
    print("PART 5: 3-CYCLE REVERSAL AS 3-ARC FLIP")
    print("  Reversing a 3-cycle = flipping 3 arcs simultaneously")
    print("  This is the 'Reidemeister R3' move")
    print("  Already verified: H-invariant at n=4, NOT at n>=5")
    print("  New question: Delta_3 on 3-cycle arc triples?")
    print("=" * 70)

    n = 5
    arcs = arcs_list(n)
    m = len(arcs)

    print(f"\n--- n = {n} ---")
    cycle_delta3 = Counter()
    noncycle_delta3 = Counter()
    total_cycle = 0
    total_noncycle = 0

    for bits in range(min(2**(n*(n-1)//2), 200)):
        A = bits_to_adj(bits, n)

        for a, b, c in combinations(range(n), 3):
            # Get the 3 arcs in this triple
            triple_arcs = [(min(a,b), max(a,b)), (min(a,c), max(a,c)), (min(b,c), max(b,c))]

            # Check if this triple forms a 3-cycle in A
            is_cycle = (A[a][b] + A[b][c] + A[c][a] == 3) or (A[b][a] + A[c][b] + A[a][c] == 3)

            delta3 = finite_diff(A, n, triple_arcs, compute_H_dp)

            if is_cycle:
                cycle_delta3[delta3] += 1
                total_cycle += 1
            else:
                noncycle_delta3[delta3] += 1
                total_noncycle += 1

    print(f"  3-cycle triples: Delta_3 distribution: {dict(sorted(cycle_delta3.items()))}")
    print(f"  Non-cycle triples: Delta_3 distribution: {dict(sorted(noncycle_delta3.items()))}")
    print(f"  3-cycle Delta_3 all zero: {all(d == 0 for d in cycle_delta3.keys())}")

    # ====================================
    # PART 6: Omega palindrome analysis
    # ====================================
    print("\n" + "=" * 70)
    print("PART 6: OMEGA DIMENSION PALINDROME (POINCARE DUALITY)")
    print("  For Paley T_7: Omega = [7,21,42,63,63,42,21] — palindromic!")
    print("  For transitive: Omega = [n, C(n,2), ..., 1] = Pascal row")
    print("  Question: which tournaments have palindromic Omega?")
    print("=" * 70)

    def compute_omega(A, n):
        paths = {0: [(v,) for v in range(n)]}
        for d in range(1, n):
            pd = []
            for p in paths[d-1]:
                last = p[-1]
                used = set(p)
                for v in range(n):
                    if v not in used and A[last][v] == 1:
                        pd.append(p + (v,))
            paths[d] = pd
            if not pd:
                break

        omega_dims = []
        allowed_set = {}
        for d in range(n):
            if d not in paths:
                omega_dims.append(0)
                allowed_set[d] = set()
                continue
            allowed_set[d] = set(paths[d])
            n_d = len(paths[d])

            if d == 0:
                omega_dims.append(n_d)
                continue

            junk_set = set()
            for p in paths[d]:
                for fi in range(d + 1):
                    face = p[:fi] + p[fi+1:]
                    if face not in allowed_set.get(d-1, set()):
                        junk_set.add(face)

            if not junk_set:
                omega_dims.append(n_d)
                continue

            junk_list = sorted(junk_set)
            C = np.zeros((len(junk_list), n_d), dtype=float)
            junk_idx = {j: i for i, j in enumerate(junk_list)}

            for j_col, p in enumerate(paths[d]):
                for fi in range(d + 1):
                    face = p[:fi] + p[fi+1:]
                    sign = 1 if fi % 2 == 0 else -1
                    if face in junk_idx:
                        C[junk_idx[face], j_col] += sign

            rank_c = int(np.linalg.matrix_rank(C)) if len(junk_list) > 0 else 0
            omega_dims.append(n_d - rank_c)

        return omega_dims

    for n in [4, 5, 6]:
        print(f"\n--- n = {n} ---")
        palindrome_count = 0
        total = 0
        omega_by_H = defaultdict(list)

        for bits in range(min(2**(n*(n-1)//2), 2000)):
            A = bits_to_adj(bits, n)
            H = compute_H_dp(A, n)
            omega = compute_omega(A, n)
            total += 1

            # Check palindrome
            is_palin = (omega == omega[::-1])
            if is_palin:
                palindrome_count += 1

            omega_by_H[H].append((tuple(omega), is_palin))

        print(f"  Palindromic Omega: {palindrome_count}/{total}")
        print(f"\n  H -> Omega palindromy:")
        for H in sorted(omega_by_H.keys()):
            items = omega_by_H[H]
            n_palin = sum(1 for _, p in items if p)
            omegas = set(o for o, _ in items)
            print(f"    H={H:3d}: {n_palin}/{len(items)} palindromic, {len(omegas)} distinct Omega")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)

if __name__ == '__main__':
    main()
