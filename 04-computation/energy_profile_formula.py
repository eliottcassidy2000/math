#!/usr/bin/env python3
"""
energy_profile_formula.py -- Finding H = f(E_2, E_3, ...) explicitly

THEOREM CANDIDATE: H(sigma) for circulant tournaments on Z_p is determined
by the representation profile r_S, equivalently by (E_2, E_3, ...).

At p=7: H = 241.5 - 3.5*E_2 (exactly linear)
At p=11: H = f(E_2, E_3) with 4 points
At p=13: H = f(E_2, E_3) with 6 points (E_2 alone fails)

This script:
1. Finds the explicit polynomial H = a + b*E_2 + c*E_3 + ... at each p
2. Tests whether it extends to p=17 (if feasible)
3. Investigates the deeper structural reason

The key insight: by OCF, H = I(Omega(T), 2) = sum_{k>=0} alpha_k * 2^k
where alpha_k = # of k-element independent sets in Omega(T).
Each alpha_k depends on the cycle structure, which for circulant T
depends on S, which groups by representation profile.

Author: kind-pasteur-2026-03-12-S60
"""

import cmath
import math
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def representation_function(S, p):
    r = [0] * p
    for a in S:
        for b in S:
            r[(a + b) % p] += 1
    return r


def energy_moments(r):
    """Compute E_k = sum r(n)^k for k=2,3,4,5."""
    return {k: sum(x**k for x in r) for k in [2, 3, 4, 5]}


def additive_energy(S, p):
    S_set = set(S)
    e = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    e += 1
    return e


def cycle_counts(A, p, max_k=None):
    """Count directed odd cycles by length."""
    if max_k is None:
        max_k = p
    counts = {}
    for k in range(3, max_k + 1, 2):
        count = 0
        for subset in combinations(range(p), k):
            verts = list(subset)
            nc = count_directed_ham_cycles_subset(A, verts)
            count += nc
        counts[k] = count
    return counts


def count_directed_ham_cycles_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        fwd = A[a][b] * A[b][c] * A[c][a]
        bwd = A[a][c] * A[c][b] * A[b][a]
        return fwd + bwd
    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def main():
    print("=" * 70)
    print("ENERGY PROFILE -> H FORMULA")
    print("=" * 70)

    # PART 1: Collect all (E_2, E_3, E_4, H) data
    print("\n--- PART 1: FULL DATA TABLE ---")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        print(f"\n  p={p}, m={m}:")
        print(f"    {'E_2':>6} {'E_3':>8} {'E_4':>10} {'H':>10} {'c3':>5} {'c5':>6}")

        # Collect one representative per (E_2, E_3) class
        class_data = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

            r = representation_function(S, p)
            EM = energy_moments(r)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)

            # Cycle counts (3-cycles always)
            c3 = 0
            for a, b, c in combinations(range(p), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3 += 1

            # 5-cycles (feasible for p <= 13)
            c5 = 0
            if p <= 13:
                for subset in combinations(range(p), 5):
                    verts = list(subset)
                    c5 += count_directed_ham_cycles_subset(A, verts)

            key = (EM[2], EM[3])
            if key not in class_data:
                class_data[key] = {'E2': EM[2], 'E3': EM[3], 'E4': EM[4], 'H': H,
                                   'c3': c3, 'c5': c5, 'S': S}

        for key in sorted(class_data):
            d = class_data[key]
            print(f"    {d['E2']:>6} {d['E3']:>8} {d['E4']:>10} {d['H']:>10} {d['c3']:>5} {d['c5']:>6}")

    # PART 2: Linear regression H = a + b*E_2 + c*E_3
    print("\n--- PART 2: H = a + b*E_2 + c*E_3 FIT ---")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        class_data = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            r = representation_function(S, p)
            EM = energy_moments(r)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            key = (EM[2], EM[3])
            if key not in class_data:
                class_data[key] = (EM[2], EM[3], EM[4], H)

        points = sorted(class_data.values())
        n = len(points)

        print(f"\n  p={p}, {n} classes:")

        # If n >= 3: solve H = a + b*E2 + c*E3 via least squares (3x3 normal equations)
        if n >= 3:
            # Construct design matrix X = [[1, E2, E3], ...]
            # and target y = [H, ...]
            X = [[1, pt[0], pt[1]] for pt in points]
            y = [pt[3] for pt in points]

            # Normal equations: X^T X beta = X^T y
            # X^T X:
            XtX = [[0]*3 for _ in range(3)]
            Xty = [0]*3
            for i in range(n):
                for j in range(3):
                    for k in range(3):
                        XtX[j][k] += X[i][j] * X[i][k]
                    Xty[j] += X[i][j] * y[i]

            # Solve 3x3 system via Cramer's rule
            def det3(M):
                return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
                       -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
                       +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

            D = det3(XtX)
            if abs(D) > 1e-10:
                # Replace columns to get a, b, c
                def replace_col(M, col, vec):
                    M2 = [row[:] for row in M]
                    for i in range(3):
                        M2[i][col] = vec[i]
                    return M2

                a = det3(replace_col(XtX, 0, Xty)) / D
                b = det3(replace_col(XtX, 1, Xty)) / D
                c = det3(replace_col(XtX, 2, Xty)) / D

                print(f"    H = {a:.4f} + {b:.4f}*E_2 + {c:.4f}*E_3")

                # Check residuals
                residuals = [y[i] - (a + b*X[i][1] + c*X[i][2]) for i in range(n)]
                max_res = max(abs(r) for r in residuals)
                print(f"    Residuals: {[f'{r:.1f}' for r in residuals]}")
                print(f"    Max residual: {max_res:.4f}")
                print(f"    EXACT linear fit? {max_res < 0.01}")
            else:
                print(f"    Singular system (D={D:.4f})")

    # PART 3: Is H a polynomial in E_2 alone (but nonlinear)?
    print("\n--- PART 3: H AS POLYNOMIAL IN E_2 ---")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        EH_map = {}
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            E = additive_energy(S, p)
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            if E not in EH_map:
                EH_map[E] = set()
            EH_map[E].add(H)

        print(f"\n  p={p}:")
        ambiguous = False
        for E in sorted(EH_map):
            H_vals = sorted(EH_map[E])
            if len(H_vals) > 1:
                ambiguous = True
            print(f"    E={E}: H={H_vals}")
        if not ambiguous:
            print(f"    E_2 alone determines H at p={p}")
        else:
            print(f"    E_2 alone does NOT determine H at p={p}")
            print(f"    Need (E_2, E_3) pair")

    # PART 4: The c3 connection
    print("\n--- PART 4: CYCLE COUNTS AS FUNCTIONS OF E ---")
    print("c3 depends only on S (not orientation). But for circulant T,")
    print("c3 = p * |zero-sum triples| / 3 which depends on S.")
    print("Question: does c3 determine E? Or vice versa?")

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        c3E_data = defaultdict(set)
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            E = additive_energy(S, p)
            A = build_adj(p, S)

            # Count 3-cycle vertex sets
            c3 = 0
            for a, b, c in combinations(range(p), 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3 += 1

            c3E_data[c3].add(E)

        print(f"\n  p={p}:")
        for c3 in sorted(c3E_data):
            E_vals = sorted(c3E_data[c3])
            print(f"    c3={c3}: E values = {E_vals}")

    # PART 5: The number-theoretic structure of E values
    print("\n--- PART 5: STRUCTURE OF E VALUES ---")
    print("E(S) for m-element subsets of Z_p obtained by orientation choices")

    for p in [7, 11, 13, 17]:
        m = (p - 1) // 2
        pairs = [(s, p - s) for s in range(1, m + 1)]

        E_vals = set()
        for bits in range(1 << m):
            sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
            S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
            E = additive_energy(S, p)
            E_vals.add(E)

        E_sorted = sorted(E_vals)
        E_diffs = [E_sorted[i+1] - E_sorted[i] for i in range(len(E_sorted)-1)]

        print(f"\n  p={p}, m={m}:")
        print(f"    E values: {E_sorted}")
        print(f"    E diffs:  {E_diffs}")
        print(f"    E range:  [{E_sorted[0]}, {E_sorted[-1]}]")
        print(f"    # distinct: {len(E_sorted)}")

    # PART 6: Why does the representation profile determine H?
    print("\n--- PART 6: STRUCTURAL EXPLANATION ---")
    print("""
    WHY does the representation profile r_S determine H?

    For a circulant tournament T_p with connection set S:
    - The adjacency matrix A has eigenvalues lambda_k = S_hat(k)
    - H = sum over Hamiltonian paths = permanent-like quantity
    - OCF: H = I(Omega(T), 2) where Omega is the conflict graph

    The conflict graph Omega depends on cycle structure.
    The cycle count c_k = (directed k-cycles) depends on the eigenvalues:
      c_k = sum_{k} lambda_{k_1} * ... * lambda_{k_k} (trace of A^k, restricted)

    Actually, for UNDIRECTED cycle vertex sets:
      c_k is determined by the Fourier magnitudes |S_hat(k)|^2 = Q_k
      (because cycle counts depend on |lambda_k|^2, not phase)

    The sorted Q_k profile is exactly the sorted |S_hat(k)|^2 profile,
    which is determined by the representation function r_S via Parseval:
      |S_hat(k)|^2 = sum_n r_S(n) * omega^{kn}  (Fourier transform of r_S)

    So: r_S determines {Q_k} determines {c_k} determines Omega determines H.

    But WAIT: two different r_S profiles can give the same SORTED Q_k set
    but H values differ because the ASSIGNMENT of Q_k to k matters for
    higher-order cycle overlap structure.

    Actually NO: the key insight is that for circulant tournaments,
    the automorphism group Z_p acts on everything, so the SORTED profile
    suffices. The overlap structure of cycles is determined by the
    MULTISET {Q_k}, not the individual assignments.

    CORRECTION: At p=13, we saw that E_2=118 has TWO different H values.
    These correspond to TWO different sorted r_S profiles with the same E_2.
    Adding E_3 distinguishes them. So it's the PROFILE that matters,
    not just a few moments.

    CONJECTURE (HYP-575): For circulant tournaments on Z_p,
    H(T) is determined by the sorted representation profile
    {r_S(n)}_{n=1}^{p-1} (multiset).
    """)

    # PART 7: Check at p=17 whether (E2, E3, E4) suffice
    print("--- PART 7: p=17 CHECK (E2, E3) ---")
    p = 17
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    print(f"  p={p}, m={m}, 2^m={1 << m} orientations (may be slow for H)")

    import time

    class_data = {}
    t0 = time.time()
    for bits in range(1 << m):
        sigma = tuple((1 if bits & (1 << i) else -1) for i in range(m))
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))

        r = representation_function(S, p)
        EM = energy_moments(r)
        key = (EM[2], EM[3])

        if key not in class_data:
            # Only compute H for one representative per class
            A = build_adj(p, S)
            H = count_ham_paths(A, p)
            class_data[key] = H
            elapsed = time.time() - t0
            print(f"    New class ({len(class_data)}): E2={EM[2]}, E3={EM[3]}, H={H} "
                  f"({elapsed:.1f}s elapsed)")

    print(f"\n  Total distinct (E2,E3) classes at p=17: {len(class_data)}")

    # Check if E2 alone suffices
    E2_groups = defaultdict(set)
    for (e2, e3), H in class_data.items():
        E2_groups[e2].add(H)

    E2_ambiguous = any(len(v) > 1 for v in E2_groups.values())
    print(f"  E2 alone determines H? {'NO' if E2_ambiguous else 'YES'}")
    if E2_ambiguous:
        for e2, H_vals in sorted(E2_groups.items()):
            if len(H_vals) > 1:
                print(f"    E2={e2}: H values = {sorted(H_vals)}")


if __name__ == '__main__':
    main()
