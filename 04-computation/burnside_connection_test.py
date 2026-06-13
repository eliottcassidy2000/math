#!/usr/bin/env python3
"""
Test whether the Even-Odd Split Lemma can be understood through
group-theoretic / Burnside-type character sums.

The alternating sum: sum_{S subset M} (-1)^|S| Delta(S, M\S) = 0
where Delta(S,R) = E_i(S)*B_j(R) - E_j(S)*B_i(R).

This looks like the sign character of Z_2^m acting on subsets of M.

Key questions:
1. Can we decompose Delta(S,R) as a character sum?
2. Is there a natural group action that explains the vanishing?
3. Does the Q=AB, K=BA factorization from Feng's paper appear?

Instance: opus-2026-03-05-S4b
"""
import random
from itertools import permutations, combinations
from math import factorial
import numpy as np

def h_end(T, n, S_list, v):
    """Count Hamiltonian paths in T[S_list union {v}] ending at v."""
    S = set(S_list)
    if not S:
        return 1
    S_sorted = sorted(S)
    m = len(S_sorted)
    idx = {s: i for i, s in enumerate(S_sorted)}
    full = (1 << m) - 1
    dp = [[0]*m for _ in range(1 << m)]
    for i, s in enumerate(S_sorted):
        dp[1 << i][i] = 1  # start at s
    for mask in range(1, (1 << m)):
        for li in range(m):
            c = dp[mask][li]
            if c == 0: continue
            for ni in range(m):
                if mask & (1 << ni): continue
                if T[S_sorted[li]*n + S_sorted[ni]]:
                    dp[mask | (1 << ni)][ni] += c
    total = 0
    for li in range(m):
        if T[S_sorted[li]*n + v]:
            total += dp[full][li]
    return total

def h_start(T, n, R_list, v):
    """Count Hamiltonian paths in T[{v} union R_list] starting at v."""
    R = set(R_list)
    if not R:
        return 1
    R_sorted = sorted(R)
    m = len(R_sorted)
    full = (1 << m) - 1
    dp = [[0]*m for _ in range(1 << m)]
    for i, r in enumerate(R_sorted):
        if T[v*n + r]:
            dp[1 << i][i] = 1
    for mask in range(1, (1 << m)):
        for li in range(m):
            c = dp[mask][li]
            if c == 0: continue
            for ni in range(m):
                if mask & (1 << ni): continue
                if T[R_sorted[li]*n + R_sorted[ni]]:
                    dp[mask | (1 << ni)][ni] += c
    return sum(dp[full])

def random_tournament(n):
    T = [0]*(n*n)
    for a in range(n):
        for b in range(a+1, n):
            if random.random() < 0.5:
                T[a*n+b] = 1
            else:
                T[b*n+a] = 1
    return T

# ==============================================================
# Test 1: Decompose the alternating sum by subset size
# ==============================================================
print("=" * 60)
print("TEST 1: Alternating sum decomposition by |S|")
print("=" * 60)

for n in [5, 6, 7]:
    print(f"\n--- n={n} ---")
    random.seed(42)

    for trial in range(3):
        T = random_tournament(n)
        i, j = 0, 1
        M = [v for v in range(n) if v != i and v != j]
        m = len(M)

        # Compute Delta(S, R) for each S
        sums_by_size = [0.0] * (m + 1)
        cross_Ei_Bj = [0.0] * (m + 1)
        cross_Ej_Bi = [0.0] * (m + 1)

        for smask in range(1 << m):
            S = [M[k] for k in range(m) if smask & (1 << k)]
            R = [M[k] for k in range(m) if not (smask & (1 << k))]
            sz = len(S)

            Ei = h_end(T, n, S, i)
            Ej = h_end(T, n, S, j)
            Bi = h_start(T, n, R, i)
            Bj = h_start(T, n, R, j)

            delta = Ei * Bj - Ej * Bi
            sign = (-1) ** sz
            sums_by_size[sz] += sign * delta
            cross_Ei_Bj[sz] += sign * Ei * Bj
            cross_Ej_Bi[sz] += sign * Ej * Bi

        total = sum(sums_by_size)
        print(f"  Trial {trial}: total alternating sum = {total}")
        print(f"    By |S|: {[int(x) for x in sums_by_size]}")
        print(f"    sum (-1)^|S| Ei*Bj = {sum(cross_Ei_Bj):.0f}")
        print(f"    sum (-1)^|S| Ej*Bi = {sum(cross_Ej_Bi):.0f}")

# ==============================================================
# Test 2: Transfer matrix / AB = BA structure
# ==============================================================
print("\n" + "=" * 60)
print("TEST 2: Transfer matrix structure (A*B vs B*A)")
print("=" * 60)
print("Define A(S) = (E_i(S), E_j(S)) and B(R) = (B_i(R), B_j(R))")
print("Delta(S,R) = det([A(S), B(R)]) = Ei*Bj - Ej*Bi")
print("")
print("Check if there's a factorization structure...")

for n in [5, 6]:
    print(f"\n--- n={n} ---")
    random.seed(42)
    T = random_tournament(n)
    i, j = 0, 1
    M = [v for v in range(n) if v != i and v != j]
    m = len(M)

    # Build the "A matrix": rows indexed by subsets S, columns = (E_i(S), E_j(S))
    # Build the "B matrix": rows indexed by subsets R, columns = (B_i(R), B_j(R))

    num_subsets = 1 << m
    A_mat = np.zeros((num_subsets, 2))  # (E_i, E_j) for each S
    B_mat = np.zeros((num_subsets, 2))  # (B_i, B_j) for each R

    for smask in range(num_subsets):
        S = [M[k] for k in range(m) if smask & (1 << k)]
        R = [M[k] for k in range(m) if not (smask & (1 << k))]

        A_mat[smask, 0] = h_end(T, n, S, i)
        A_mat[smask, 1] = h_end(T, n, S, j)
        # B indexed by complement
        cmask = ((1 << m) - 1) ^ smask
        B_mat[cmask, 0] = h_start(T, n, R, i)
        B_mat[cmask, 1] = h_start(T, n, R, j)

    # The alternating sum is: sum_S (-1)^|S| (A[S,0]*B[S^c,1] - A[S,1]*B[S^c,0])
    # = sum_S (-1)^|S| det([A[S], B[S^c]])

    # Sign vector
    signs = np.array([(-1)**bin(s).count('1') for s in range(num_subsets)])

    # Compute A^T * diag(signs) * B_comp
    # where B_comp[s] = B[s^c]
    B_comp = np.zeros_like(B_mat)
    full_mask = (1 << m) - 1
    for s in range(num_subsets):
        B_comp[s] = B_mat[full_mask ^ s]

    # Matrix product: A^T * D * B_comp where D = diag(signs)
    M_transfer = A_mat.T @ np.diag(signs) @ B_comp
    print(f"  A^T * diag(signs) * B_comp =")
    print(f"    [[{M_transfer[0,0]:.0f}, {M_transfer[0,1]:.0f}],")
    print(f"     [{M_transfer[1,0]:.0f}, {M_transfer[1,1]:.0f}]]")
    print(f"  (The alternating sum = M[0,1] - M[1,0] = {M_transfer[0,1] - M_transfer[1,0]:.0f})")

    # Check: is M_transfer symmetric? (would explain vanishing of off-diagonal diff)
    print(f"  M symmetric? M[0,1]={M_transfer[0,1]:.0f}, M[1,0]={M_transfer[1,0]:.0f}")
    print(f"  -> {'YES' if abs(M_transfer[0,1] - M_transfer[1,0]) < 0.5 else 'NO'} (diff={M_transfer[0,1]-M_transfer[1,0]:.0f})")

    # Also check diagonal
    print(f"  M[0,0]={M_transfer[0,0]:.0f}, M[1,1]={M_transfer[1,1]:.0f}")

# ==============================================================
# Test 3: Character decomposition
# ==============================================================
print("\n" + "=" * 60)
print("TEST 3: Fourier/character decomposition of Delta(S,R)")
print("=" * 60)
print("Expand Delta(S,R) in the Fourier basis of Z_2^m")
print("chi_T(S) = (-1)^|S cap T| for subset T of M")

for n in [5, 6]:
    print(f"\n--- n={n} ---")
    random.seed(42)
    T_tour = random_tournament(n)
    i, j = 0, 1
    M = [v for v in range(n) if v != i and v != j]
    m = len(M)

    # Compute Delta(S, M\S) for all S
    delta_vals = np.zeros(1 << m)
    for smask in range(1 << m):
        S = [M[k] for k in range(m) if smask & (1 << k)]
        R = [M[k] for k in range(m) if not (smask & (1 << k))]
        Ei = h_end(T_tour, n, S, i)
        Ej = h_end(T_tour, n, S, j)
        Bi = h_start(T_tour, n, R, i)
        Bj = h_start(T_tour, n, R, j)
        delta_vals[smask] = Ei * Bj - Ej * Bi

    # Fourier transform: f_hat(T) = (1/2^m) sum_S (-1)^|S cap T| f(S)
    fourier = np.zeros(1 << m)
    for tmask in range(1 << m):
        s = 0.0
        for smask in range(1 << m):
            inter = bin(smask & tmask).count('1')
            s += ((-1)**inter) * delta_vals[smask]
        fourier[tmask] = s / (1 << m)

    # The sign character corresponds to tmask = full_mask (all bits set)
    # (-1)^|S cap M| = (-1)^|S|
    full = (1 << m) - 1
    print(f"  Fourier coeff at full mask (sign char): {fourier[full]:.4f}")
    print(f"  Fourier coeff at 0 (trivial char): {fourier[0]:.4f}")

    # Find all nonzero Fourier coefficients
    nonzero = [(tmask, fourier[tmask]) for tmask in range(1 << m)
               if abs(fourier[tmask]) > 0.01]
    print(f"  Nonzero Fourier coeffs: {len(nonzero)} / {1 << m}")
    if len(nonzero) <= 10:
        for tmask, val in nonzero:
            bits = bin(tmask).count('1')
            print(f"    mask={tmask:0{m}b} (weight {bits}): {val:.2f}")

    # KEY: the sign character coefficient is the alternating sum / 2^m
    # Our lemma says this is 0. What about other characters?
    # Group by Hamming weight (= character of Z_2^m subgroup)
    by_weight = {}
    for tmask in range(1 << m):
        w = bin(tmask).count('1')
        by_weight.setdefault(w, 0)
        by_weight[w] += fourier[tmask]
    print(f"  Sum of Fourier coeffs by weight: {dict(sorted(by_weight.items()))}")
