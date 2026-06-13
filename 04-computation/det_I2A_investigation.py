#!/usr/bin/env python3
"""
det_I2A_investigation.py -- kind-pasteur-2026-03-13-S61

REMARKABLE FINDING: For the two deepest-ambiguity tournaments:
  bits=4728 (H=109, c7=8): det(I+2A) = 169 = 13^2
  bits=4658 (H=111, c7=9): det(I+2A) = 361 = 19^2

Also: det(I+A) = 8 vs 10 (difference = 2 = delta_c7)

QUESTIONS:
1. Is det(I+2A) always a perfect square for tournaments?
2. What is the relationship between det(I+2A) and H?
3. Does det(I+A) = c7_dir + constant?
4. Is det(I+xA) related to the independence polynomial I(Omega,x)?
5. What do the Pfaffian principal minors {49,1,9,49,1,9,1} vs {81,1,9,9,81,1,1} tell us?

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k <= 2:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def det_small(M, n):
    """Determinant via cofactor expansion for n <= 8."""
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]
    result = 0
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [[M[r][c] for c in range(n) if c != j] for r in range(1, n)]
        cofactor = ((-1) ** j) * M[0][j] * det_small(minor, n-1)
        result += cofactor
    return result


def char_poly_eval(A, n, x_val):
    """Evaluate det(xI - A) at a specific x value."""
    M = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            M[i][j] = -A[i][j]
        M[i][i] = x_val - A[i][i]
    return det_small(M, n)


# ========================================================================
# ANALYSIS 1: det(I+xA) for various x and both target tournaments
# ========================================================================
n = 7
print("=" * 70)
print("ANALYSIS 1: det(I+xA) AS A FUNCTION OF x")
print("=" * 70)

targets = [(4728, "H=109, c7=8"), (4658, "H=111, c7=9")]

for bits, label in targets:
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    c7 = count_directed_ham_cycles_on_subset(A, list(range(n)))

    print(f"\n  {label}:")
    print(f"    H={H}, c7={c7}")

    for x in range(-3, 8):
        M = [[(x * A[i][j]) for j in range(n)] for i in range(n)]
        for i in range(n):
            M[i][i] += 1
        d = det_small(M, n)
        # Check if perfect square
        is_sq = ""
        ad = abs(d)
        r = int(ad**0.5 + 0.5)
        if r*r == ad:
            is_sq = f" = {'-' if d < 0 else ''}{r}^2"
        print(f"    det(I+{x}A) = {d}{is_sq}")


# ========================================================================
# ANALYSIS 2: det(I+A) and det(I+2A) across ALL n=5 tournaments
# ========================================================================
print(f"\n\n{'='*70}")
print("ANALYSIS 2: det(I+A) AND det(I+2A) AT n=5 (EXHAUSTIVE)")
print("=" * 70)

n5 = 5
m5 = n5 * (n5 - 1) // 2
total5 = 1 << m5

by_H_n5 = defaultdict(list)
for bits in range(total5):
    A = binary_to_tournament(bits, n5)
    H = count_ham_paths(A, n5)

    M1 = [[A[i][j] for j in range(n5)] for i in range(n5)]
    for i in range(n5):
        M1[i][i] += 1
    d1 = det_small(M1, n5)

    M2 = [[2*A[i][j] for j in range(n5)] for i in range(n5)]
    for i in range(n5):
        M2[i][i] += 1
    d2 = det_small(M2, n5)

    c5 = count_directed_ham_cycles_on_subset(A, list(range(n5)))

    by_H_n5[H].append((bits, d1, d2, c5))

print(f"  H values: {sorted(by_H_n5.keys())}")
for H_val in sorted(by_H_n5.keys()):
    group = by_H_n5[H_val]
    d1_vals = sorted(set(d1 for _, d1, _, _ in group))
    d2_vals = sorted(set(d2 for _, _, d2, _ in group))
    c5_vals = sorted(set(c5 for _, _, _, c5 in group))
    print(f"  H={H_val}: det(I+A)={d1_vals}, det(I+2A)={d2_vals}, c5={c5_vals}")


# ========================================================================
# ANALYSIS 3: Does det(I+2A) relate to I(Omega, 2) = H?
# ========================================================================
print(f"\n\n{'='*70}")
print("ANALYSIS 3: RELATIONSHIP det(I+2A) vs H AT n=5")
print("=" * 70)

for H_val in sorted(by_H_n5.keys()):
    group = by_H_n5[H_val]
    for bits, d1, d2, c5 in group[:3]:  # sample
        print(f"  H={H_val}, bits={bits}: det(I+A)={d1}, det(I+2A)={d2}, c5={c5}, d2-H={d2-H_val}")


# ========================================================================
# ANALYSIS 4: Exhaustive at n=6
# ========================================================================
print(f"\n\n{'='*70}")
print("ANALYSIS 4: det(I+A) AND det(I+2A) AT n=6 (EXHAUSTIVE)")
print("=" * 70)

n6 = 6
m6 = n6 * (n6 - 1) // 2
total6 = 1 << m6

by_H_n6 = defaultdict(list)
for bits in range(total6):
    A = binary_to_tournament(bits, n6)
    H = count_ham_paths(A, n6)

    M1 = [[A[i][j] for j in range(n6)] for i in range(n6)]
    for i in range(n6):
        M1[i][i] += 1
    d1 = det_small(M1, n6)

    M2 = [[2*A[i][j] for j in range(n6)] for i in range(n6)]
    for i in range(n6):
        M2[i][i] += 1
    d2 = det_small(M2, n6)

    by_H_n6[H].append((bits, d1, d2))

print(f"  H values: {sorted(by_H_n6.keys())}")

# Check: does det(I+2A) determine H?
detH_map = defaultdict(set)
for H_val, group in by_H_n6.items():
    for bits, d1, d2 in group:
        detH_map[d2].add(H_val)

ambig = sum(1 for d2, Hs in detH_map.items() if len(Hs) > 1)
print(f"  det(I+2A) -> H ambiguities: {ambig} out of {len(detH_map)} det values")

# Show the ambiguous cases
if ambig > 0:
    for d2, Hs in sorted(detH_map.items()):
        if len(Hs) > 1:
            print(f"    det(I+2A)={d2}: H={sorted(Hs)}")


# ========================================================================
# ANALYSIS 5: det(I+2A) at n=7 for the lambda-ambiguous class
# ========================================================================
print(f"\n\n{'='*70}")
print("ANALYSIS 5: det(I+2A) AT n=7 FOR LAMBDA-AMBIGUOUS CLASS")
print("=" * 70)

# We already know bits 4728 and 4658 are in the ambiguous class
# Check: does det(I+2A) distinguish them?
for bits, label in targets:
    A = binary_to_tournament(bits, n)
    M2 = [[2*A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        M2[i][i] += 1
    d2 = det_small(M2, n)

    # Also compute det(I-A) and det(I-2A)
    Mm1 = [[-A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        Mm1[i][i] += 1
    dm1 = det_small(Mm1, n)

    Mm2 = [[-2*A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        Mm2[i][i] += 1
    dm2 = det_small(Mm2, n)

    print(f"  {label}:")
    print(f"    det(I+2A) = {d2}")
    print(f"    det(I-A) = {dm1}")
    print(f"    det(I-2A) = {dm2}")

    # The skew-adjacency matrix S = A - A^T
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

    # det(I+xS) for various x
    for x in [1, 2]:
        M = [[(x * S[i][j]) for j in range(n)] for i in range(n)]
        for i in range(n):
            M[i][i] += 1
        d = det_small(M, n)
        print(f"    det(I+{x}S) = {d}")

    # Characteristic polynomial of A at key points
    for x_val in [0, 1, 2, -1, -2]:
        cp = char_poly_eval(A, n, x_val)
        print(f"    char_poly(A, {x_val}) = det({x_val}I - A) = {cp}")


# ========================================================================
# ANALYSIS 6: Is there a formula det(I+2A) = f(H)?
# ========================================================================
print(f"\n\n{'='*70}")
print("ANALYSIS 6: FORMULA SEARCH det(I+2A) vs H AT n=5")
print("=" * 70)

print(f"\n  At n=5:")
for H_val in sorted(by_H_n5.keys()):
    group = by_H_n5[H_val]
    d2_vals = sorted(set(d2 for _, _, d2, _ in group))
    d1_vals = sorted(set(d1 for _, d1, _, _ in group))
    c5_vals = sorted(set(c5 for _, _, _, c5 in group))
    # Check d2 = aH + b?
    if len(d2_vals) == 1:
        print(f"  H={H_val}: det(I+2A)={d2_vals[0]}, det(I+A)={d1_vals[0]}, c5={c5_vals}")
    else:
        print(f"  H={H_val}: det(I+2A)={d2_vals} (MULTIPLE), det(I+A)={d1_vals}")


# ========================================================================
# ANALYSIS 7: The Pfaffian structure
# ========================================================================
print(f"\n\n{'='*70}")
print("ANALYSIS 7: PFAFFIAN PRINCIPAL MINORS")
print("=" * 70)

# For the skew-adjacency S = A - A^T (antisymmetric),
# det(S[V\v]) = Pf(S[V\v])^2 for even-size matrices.
# At n=7 (odd), the full det(S) = 0, but 6x6 minors are nonzero.

# The multiset of Pfaffian squares characterizes something...
for bits, label in targets:
    A = binary_to_tournament(bits, n)
    S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

    print(f"\n  {label}:")

    pf_squares = []
    for del_v in range(n):
        remaining = [i for i in range(n) if i != del_v]
        S6 = [[S[remaining[i]][remaining[j]] for j in range(6)] for i in range(6)]
        det6 = det_small(S6, 6)
        pf_squares.append(det6)

    print(f"    Pfaffian^2 per deleted vertex: {pf_squares}")
    print(f"    Sorted: {sorted(pf_squares)}")
    print(f"    Sum: {sum(pf_squares)}")
    print(f"    Product: {1}")  # placeholder
    for i, p in enumerate(pf_squares):
        r = int(abs(p)**0.5 + 0.5)
        print(f"      del v={i}: Pf^2={p}, Pf={r}")

    # The Pfaffian of the skew-adjacency matrix counts...
    # For a tournament, Pf(S[V\v]) counts the number of
    # "near-perfect matchings" avoiding vertex v, with signs.
    # This is related to the cofactor expansion of the determinant.

    # Also: for skew-symmetric matrices, det = Pf^2
    # and Pf counts "directed perfect matchings" (dimers) with signs.

    # Check: sum of Pfaffians = ?
    pfs = [int(abs(p)**0.5 + 0.5) for p in pf_squares]
    print(f"    Pfaffian values: {pfs}")
    print(f"    Sum of Pfaffians: {sum(pfs)}")

# ========================================================================
# ANALYSIS 8: det(I+2A) and permanent
# ========================================================================
print(f"\n\n{'='*70}")
print("ANALYSIS 8: PERMANENT OF A")
print("=" * 70)

def permanent_small(M, n):
    """Permanent via Ryser's formula."""
    total = 0
    for S in range(1, 1 << n):
        sign = (-1) ** (n - bin(S).count('1'))
        col_sums_prod = 1
        for j in range(n):
            s = 0
            for i in range(n):
                if S & (1 << i):
                    s += M[i][j]
            col_sums_prod *= s
        total += sign * col_sums_prod
    return total * ((-1) ** n)


for bits, label in targets:
    A = binary_to_tournament(bits, n)
    perm_A = permanent_small(A, n)

    # Also perm(I+A)
    IA = [[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        IA[i][i] += 1
    perm_IA = permanent_small(IA, n)

    print(f"\n  {label}:")
    print(f"    perm(A) = {perm_A}")
    print(f"    perm(I+A) = {perm_IA}")
    print(f"    H = {count_ham_paths(A, n)}")
    print(f"    perm(I+A) - H = {perm_IA - count_ham_paths(A, n)}")


print(f"\n\n{'='*70}")
print("DONE.")
print("=" * 70)
