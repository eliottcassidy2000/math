#!/usr/bin/env python3
"""
CRITICAL TEST: Is the transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) SYMMETRIC?

If M[i][j] = M[j][i] always, then:
  sum (-1)^|S| E_i(S)*B_j(R) = sum (-1)^|S| E_j(S)*B_i(R)
and the Even-Odd Split Lemma (their difference = 0) follows immediately.

The symmetry M = M^T is a STRONGER statement than just the off-diagonal difference vanishing.
It also constrains the diagonal: M[i][i] and M[j][j] independently.

Instance: opus-2026-03-05-S4b
"""
import random

def h_end(T, n, S_list, v):
    S = set(S_list)
    if not S:
        return 1
    S_sorted = sorted(S)
    m = len(S_sorted)
    full = (1 << m) - 1
    dp = [[0]*m for _ in range(1 << m)]
    for k, s in enumerate(S_sorted):
        dp[1 << k][k] = 1
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
    R = set(R_list)
    if not R:
        return 1
    R_sorted = sorted(R)
    m = len(R_sorted)
    full = (1 << m) - 1
    dp = [[0]*m for _ in range(1 << m)]
    for k, r in enumerate(R_sorted):
        if T[v*n + r]:
            dp[1 << k][k] = 1
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

def compute_transfer_matrix(T, n, i, j):
    M_verts = [v for v in range(n) if v != i and v != j]
    m = len(M_verts)
    result = [[0, 0], [0, 0]]
    for smask in range(1 << m):
        S = [M_verts[k] for k in range(m) if smask & (1 << k)]
        R = [M_verts[k] for k in range(m) if not (smask & (1 << k))]
        sz = len(S)
        sign = (-1) ** sz
        Ei = h_end(T, n, S, i)
        Ej = h_end(T, n, S, j)
        Bi = h_start(T, n, R, i)
        Bj = h_start(T, n, R, j)
        result[0][0] += sign * Ei * Bi
        result[0][1] += sign * Ei * Bj
        result[1][0] += sign * Ej * Bi
        result[1][1] += sign * Ej * Bj
    return result

print("SYMMETRY CHECK: Is M[i][j] = M[j][i] always?")
print("=" * 60)

for n in range(4, 9):
    random.seed(42)
    trials = {4: 500, 5: 500, 6: 200, 7: 50, 8: 10}[n]
    sym_count = 0
    total_tests = 0
    fails = []

    for trial in range(trials):
        T = random_tournament(n)
        for i in range(min(n, 3)):
            for j in range(i+1, min(n, 4)):
                M = compute_transfer_matrix(T, n, i, j)
                total_tests += 1
                if M[0][1] == M[1][0]:
                    sym_count += 1
                else:
                    if len(fails) < 3:
                        fails.append((trial, i, j, M))

    if sym_count == total_tests:
        print(f"n={n}: {sym_count}/{total_tests} SYMMETRIC (100%)")
    else:
        print(f"n={n}: {sym_count}/{total_tests} symmetric ({100*sym_count/total_tests:.1f}%)")
        for trial, i, j, M in fails:
            print(f"  FAIL trial={trial} ({i},{j}): M={M}")

# Additional: check diagonal pattern
print("\n" + "=" * 60)
print("DIAGONAL ANALYSIS: What is M[i][i] + M[j][j]?")
print("=" * 60)

for n in [5, 6, 7]:
    random.seed(42)
    trials = 100 if n <= 6 else 30
    diag_sums = set()

    for trial in range(trials):
        T = random_tournament(n)
        i, j = 0, 1
        M = compute_transfer_matrix(T, n, i, j)
        dsum = M[0][0] + M[1][1]
        diag_sums.add(dsum)

    print(f"n={n}: M[i][i]+M[j][j] values seen: {sorted(diag_sums)}")

# Check M[i][i] - M[j][j]
print("\nM[i][i] - M[j][j] values:")
for n in [5, 6, 7]:
    random.seed(42)
    trials = 100 if n <= 6 else 30
    diag_diffs = {}

    for trial in range(trials):
        T = random_tournament(n)
        i, j = 0, 1
        M = compute_transfer_matrix(T, n, i, j)
        ddiff = M[0][0] - M[1][1]
        diag_diffs[ddiff] = diag_diffs.get(ddiff, 0) + 1

    print(f"n={n}: {dict(sorted(diag_diffs.items()))}")

# Check: M[i][i] + M[j][j] = ?  Maybe related to H(T)?
print("\n" + "=" * 60)
print("RELATING DIAGONAL TO H(T)")
print("=" * 60)

def ham_count(T, n):
    FULL = (1 << n) - 1
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for last in range(n):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(n):
                if mask & (1<<nxt): continue
                if T[last*n+nxt]:
                    dp[mask|(1<<nxt)][nxt] += c
    return sum(dp[FULL])

for n in [5, 6, 7]:
    random.seed(42)
    trials = 50 if n <= 6 else 20
    print(f"\nn={n}:")
    for trial in range(min(5, trials)):
        T = random_tournament(n)
        i, j = 0, 1
        M = compute_transfer_matrix(T, n, i, j)
        H = ham_count(T, n)
        adj_ij = sum(1 for P_unused in [None])  # placeholder
        # Count paths with i before j, and j before i
        # M[i][i] = sum (-1)^|S| E_i(S)*B_i(R)
        # M[j][j] = sum (-1)^|S| E_j(S)*B_j(R)
        trace = M[0][0] + M[1][1]
        diff = M[0][0] - M[1][1]
        offdiag = M[0][1]  # = M[1][0] by symmetry
        print(f"  H={H}, trace={trace}, diff={diff}, offdiag={offdiag}, det={M[0][0]*M[1][1]-offdiag**2}")
