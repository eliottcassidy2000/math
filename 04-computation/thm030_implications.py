#!/usr/bin/env python3
"""
IMPLICATIONS OF THM-030 (Transfer Matrix Symmetry).

Now that M[a,b] = M[b,a] is PROVED, what follows?

1. Off-diagonal sum: sum_{a!=b} M[a,b] = ? (proved here)
2. Row sum = column sum (immediate from symmetry)
3. M as a symmetric bilinear form: what does it measure?
4. Spectral properties of M (real eigenvalues, etc.)
5. Connection to H(T) via trace formula (THM-027)

kind-pasteur-2026-03-06-S25
"""

from itertools import permutations
from sympy import symbols, expand, Poly, Matrix, factor
import numpy as np
from collections import defaultdict

def setup(n):
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')
    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]
    def t(i, j):
        if i == j: return 0
        return r + s(i, j)
    return r, sv, s, t

def transfer_M(t_fn, vertex_set, a, b):
    V = sorted(vertex_set)
    U = [v for v in V if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        S_set = set(S) | {a}
        R_set = set(R) | {b}
        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1): prod *= t_fn(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1
        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1): prod *= t_fn(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1
        result += sign * ea * bb
    return expand(result)

def total_ham(vertex_set, t_fn):
    vs = sorted(vertex_set)
    if len(vs) <= 1: return 1
    total = 0
    for perm in permutations(vs):
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)

def ham_paths_ab(vertex_set, a, b, t_fn):
    vs = sorted(vertex_set)
    total = 0
    for perm in permutations(vs):
        if perm[0] != a or perm[-1] != b: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)


# ============================================================
# 1. Off-diagonal sum at r=c/2 (numeric tournaments)
# ============================================================
print("=" * 70)
print("1. Off-diagonal sum and trace of M for numeric tournaments")
print("=" * 70)

import random
random.seed(42)

for n in range(3, 8):
    # Random tournament
    T = {}
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[(i,j)] = 1; T[(j,i)] = 0
            else:
                T[(i,j)] = 0; T[(j,i)] = 1

    def t_num(i, j):
        if i == j: return 0
        return T.get((i,j), 0)

    W = set(range(n))
    trace = 0
    off_diag = 0
    M_mat = {}
    for a in range(n):
        for b in range(n):
            if a == b:
                # Compute M[a,a] (diagonal entry)
                # Use the trace formula: M[a,a] = (-1)^{pos(a,P)} summed over Ham paths
                M_aa = 0
                for perm in permutations(range(n)):
                    prod = 1
                    for k in range(len(perm)-1):
                        prod *= t_num(perm[k], perm[k+1])
                    if prod != 0:
                        # Find position of a in perm
                        pos = list(perm).index(a)
                        M_aa += (-1)**pos * prod
                M_mat[(a,a)] = M_aa
                trace += M_aa
            else:
                M_ab = 0
                U = [v for v in range(n) if v != a and v != b]
                for mask in range(1 << len(U)):
                    S = [U[k] for k in range(len(U)) if mask & (1 << k)]
                    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                    sign = (-1)**len(S)
                    S_set = set(S) | {a}
                    R_set = set(R) | {b}
                    ea = 0
                    for p in permutations(sorted(S_set)):
                        if p[-1] != a: continue
                        prod = 1
                        for k in range(len(p)-1):
                            prod *= t_num(p[k], p[k+1])
                        ea += prod
                    if len(S_set) == 1: ea = 1
                    bb = 0
                    for p in permutations(sorted(R_set)):
                        if p[0] != b: continue
                        prod = 1
                        for k in range(len(p)-1):
                            prod *= t_num(p[k], p[k+1])
                        bb += prod
                    if len(R_set) == 1: bb = 1
                    M_ab += sign * ea * bb
                M_mat[(a,b)] = M_ab
                off_diag += M_ab

    # Compute H(T) directly
    H = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= t_num(perm[k], perm[k+1])
        H += prod

    # Check symmetry
    sym = all(M_mat[(a,b)] == M_mat[(b,a)] for a in range(n) for b in range(n) if a != b)

    print(f"\n  n={n}: H(T) = {H}")
    print(f"    trace(M) = {trace}")
    print(f"    off_diag_sum = {off_diag}")
    print(f"    M symmetric: {sym}")

    if n % 2 == 1:
        print(f"    n odd: trace = {trace}, expected H = {H}, match = {trace == H}")
        print(f"    n odd: off_diag = {off_diag}, expected 0, match = {off_diag == 0}")
    else:
        print(f"    n even: trace = {trace}, expected 0, match = {trace == 0}")
        print(f"    n even: off_diag = {off_diag}, expected 2H = {2*H}, match = {off_diag == 2*H}")


# ============================================================
# 2. Eigenvalue structure of M (numeric)
# ============================================================
print("\n" + "=" * 70)
print("2. Eigenvalue structure of M")
print("=" * 70)

for n in range(3, 8):
    T = {}
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[(i,j)] = 1; T[(j,i)] = 0
            else:
                T[(i,j)] = 0; T[(j,i)] = 1

    def t_num2(i, j):
        if i == j: return 0
        return T.get((i,j), 0)

    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            if a == b:
                val = 0
                for perm in permutations(range(n)):
                    prod = 1
                    for k in range(len(perm)-1):
                        prod *= t_num2(perm[k], perm[k+1])
                    if prod != 0:
                        pos = list(perm).index(a)
                        val += (-1)**pos * prod
                M[a,a] = val
            else:
                val = 0
                U = [v for v in range(n) if v != a and v != b]
                for mask in range(1 << len(U)):
                    S = [U[k] for k in range(len(U)) if mask & (1 << k)]
                    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                    sign = (-1)**len(S)
                    S_set = set(S) | {a}
                    R_set = set(R) | {b}
                    ea = 0
                    for p in permutations(sorted(S_set)):
                        if p[-1] != a: continue
                        prod = 1
                        for k in range(len(p)-1):
                            prod *= t_num2(p[k], p[k+1])
                        ea += prod
                    if len(S_set) == 1: ea = 1
                    bb = 0
                    for p in permutations(sorted(R_set)):
                        if p[0] != b: continue
                        prod = 1
                        for k in range(len(p)-1):
                            prod *= t_num2(p[k], p[k+1])
                        bb += prod
                    if len(R_set) == 1: bb = 1
                    val += sign * ea * bb
                M[a,b] = val

    eigenvalues = np.linalg.eigvalsh(M)  # symmetric -> real eigenvalues
    H = sum(1 for perm in permutations(range(n))
            if all(T.get((perm[k], perm[k+1]), 0) == 1 for k in range(len(perm)-1)))

    print(f"\n  n={n}: H(T)={H}")
    print(f"    eigenvalues = {sorted(eigenvalues)[::-1]}")
    print(f"    sum(evals) = {sum(eigenvalues):.1f} (= trace = {'H' if n%2==1 else '0'})")
    print(f"    max eval = {max(eigenvalues):.3f}")

    # Check: is max eigenvalue related to H?
    # For symmetric M, the spectral norm is max |eigenvalue|
    print(f"    spectral norm = {max(abs(eigenvalues)):.3f}")

    # Is H/n an eigenvalue? (From trace = H for odd n, suggests average eigenvalue = H/n)
    if n % 2 == 1:
        # Check if H/n is close to an eigenvalue
        closest = min(eigenvalues, key=lambda x: abs(x - H/n))
        print(f"    H/n = {H/n:.3f}, closest eigenvalue = {closest:.3f}")


# ============================================================
# 3. M as quadratic form: v^T M v for characteristic vectors
# ============================================================
print("\n" + "=" * 70)
print("3. M applied to special vectors")
print("=" * 70)

# For the all-ones vector: sum_ab M[a,b] = trace + off_diag
# For odd n: H + 0 = H. For even n: 0 + 2H = 2H.
# So 1^T M 1 = H (odd) or 2H (even).
# More precisely: 1^T M 1 = sum_{a,b} M[a,b] = [1+(-1)^n] H

# But H = sum over Ham paths of product. Does v^T M v for other v give nice things?

# Score vector: s_i = out-degree of i
# Check: does s^T M s relate to any tournament invariant?

for n in [5, 7]:
    T = {}
    # Use the Paley-like tournament for cleaner structure
    # T_n = tournament on Z_n where i->j iff (j-i) is a QR mod n
    for i in range(n):
        for j in range(i+1, n):
            diff = (j - i) % n
            # QR mod n
            is_qr = any((k*k) % n == diff for k in range(1, n))
            if is_qr:
                T[(i,j)] = 1; T[(j,i)] = 0
            else:
                T[(i,j)] = 0; T[(j,i)] = 1

    def t_num3(i, j):
        if i == j: return 0
        return T.get((i,j), 0)

    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            if a == b:
                val = 0
                for perm in permutations(range(n)):
                    prod = 1
                    for k in range(len(perm)-1):
                        prod *= t_num3(perm[k], perm[k+1])
                    if prod != 0:
                        pos = list(perm).index(a)
                        val += (-1)**pos * prod
                M[a,a] = val
            else:
                val = 0
                U = [v for v in range(n) if v != a and v != b]
                for mask in range(1 << len(U)):
                    S = [U[k] for k in range(len(U)) if mask & (1 << k)]
                    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                    sign = (-1)**len(S)
                    S_set = set(S) | {a}
                    R_set = set(R) | {b}
                    ea = 0
                    for p in permutations(sorted(S_set)):
                        if p[-1] != a: continue
                        prod = 1
                        for k in range(len(p)-1):
                            prod *= t_num3(p[k], p[k+1])
                        ea += prod
                    if len(S_set) == 1: ea = 1
                    bb = 0
                    for p in permutations(sorted(R_set)):
                        if p[0] != b: continue
                        prod = 1
                        for k in range(len(p)-1):
                            prod *= t_num3(p[k], p[k+1])
                        bb += prod
                    if len(R_set) == 1: bb = 1
                    val += sign * ea * bb
                M[a,b] = val

    H = sum(1 for perm in permutations(range(n))
            if all(T.get((perm[k], perm[k+1]), 0) == 1 for k in range(len(perm)-1)))

    scores = np.array([sum(T.get((i,j),0) for j in range(n) if j != i) for i in range(n)])
    ones = np.ones(n)

    print(f"\n  n={n} (Paley-like): H(T) = {H}")
    print(f"    scores = {scores}")
    print(f"    M = \n{M}")
    print(f"    eigenvalues = {sorted(np.linalg.eigvalsh(M))[::-1]}")
    print(f"    1^T M 1 = {ones @ M @ ones:.1f}")
    print(f"    s^T M s = {scores @ M @ scores:.1f}")

    # For vertex-transitive tournaments: all M[a,a] equal, all M[a,b] (a!=b) equal
    diag_vals = set(int(round(M[a,a])) for a in range(n))
    off_diag_vals = set(int(round(M[a,b])) for a in range(n) for b in range(n) if a != b)
    print(f"    diagonal values: {diag_vals}")
    print(f"    off-diagonal values: {off_diag_vals}")
    if len(diag_vals) == 1 and len(off_diag_vals) == 1:
        d = list(diag_vals)[0]
        o = list(off_diag_vals)[0]
        print(f"    M = {d}*I + {o}*(J-I) = {o}*J + {d-o}*I")
        print(f"    eigenvalues: {d-o} (mult n-1), {d+(n-1)*o} (mult 1)")
        print(f"    {d-o} = {d-o}, {d+(n-1)*o} = {d+(n-1)*o}")


# ============================================================
# 4. Does symmetry + trace = H characterize M uniquely?
# ============================================================
print("\n" + "=" * 70)
print("4. Row sums of M")
print("=" * 70)

# Now that M is symmetric, row sum = column sum.
# Column sum cs(b) = sum_{a != b} M[a,b].
# For odd n: sum_b cs(b) = off_diag = 0.
# For even n: sum_b cs(b) = off_diag = 2H.

# What IS the row sum? Let's check if it equals B_b(W) or some function of it.

for n in [3, 4, 5]:
    r_sym, sv, s, t = setup(n)
    W = set(range(n))
    print(f"\n  n={n}:")
    for b in sorted(W):
        rs = 0
        for a in sorted(W - {b}):
            rs = expand(rs + transfer_M(t, W, a, b))
        # Also compute B_b
        Bb = 0
        for perm in permutations(sorted(W)):
            if perm[0] != b: continue
            prod = 1
            for k in range(len(perm)-1):
                prod *= t(perm[k], perm[k+1])
            Bb += prod
        Bb = expand(Bb)
        print(f"    b={b}: row_sum = {rs}")
        # Is row_sum even in r?
        if rs != 0:
            p = Poly(rs, r_sym)
            parities = set(d[0] % 2 for d in p.as_dict().keys())
            print(f"           r-parities: {parities}")


# ============================================================
# 5. Symbolic eigenvalue structure
# ============================================================
print("\n" + "=" * 70)
print("5. Symbolic M at n=3 (eigenvalues)")
print("=" * 70)

r_sym, sv, s, t = setup(3)
W = {0, 1, 2}
M_sym = {}
for a in sorted(W):
    for b in sorted(W):
        if a == b:
            # Diagonal: sum over paths through, with (-1)^pos
            val = 0
            for perm in permutations(sorted(W)):
                prod = 1
                for k in range(len(perm)-1):
                    prod *= t(perm[k], perm[k+1])
                pos = list(perm).index(a)
                val += (-1)**pos * prod
            M_sym[(a,a)] = expand(val)
        else:
            M_sym[(a,b)] = transfer_M(t, W, a, b)

print("  M matrix:")
for a in sorted(W):
    row = [str(M_sym[(a,b)]) for b in sorted(W)]
    print(f"    [{', '.join(row)}]")

# Characteristic polynomial
M_mat = Matrix(3, 3, lambda i, j: M_sym[(i,j)])
char_poly = M_mat.charpoly()
print(f"\n  Characteristic polynomial: {char_poly}")

# For n=3 with c-tournament specialization (r=1/2):
# M entries are s-values. char poly should be nice.

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
