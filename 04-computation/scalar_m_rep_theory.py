#!/usr/bin/env python3
"""
REPRESENTATION THEORY APPROACH to proving M = (H/n)*I for VT tournaments.

For a VT tournament T on n vertices with automorphism group G:
  - G acts transitively on [n]
  - M is G-equivariant: M[g(a), g(b)] = M[a,b] for all g in G
  - M is symmetric (THM-030): M[a,b] = M[b,a]

By Schur's lemma, if the permutation representation of G on C^n decomposes as
  C^n = V_0 + V_1 + ... + V_k
with V_0 = trivial rep = span(1,...,1), then M acts as a scalar on each V_i.

M is scalar (= lambda*I) iff ALL these scalars are equal.

For G = Z/nZ (circulant): C^n decomposes into n 1-dim irreps chi_j.
  M acts as lambda_j = sum_{k=0}^{n-1} M[0,k] omega^{jk} on chi_j.
  All lambda_j = H/n iff M[0,k] = 0 for k >= 1.

KEY QUESTION: Why does the transfer matrix M respect MORE symmetry
than just the automorphism group? What additional constraint forces
the Fourier coefficients to vanish?

APPROACH: The transfer matrix satisfies M[a,b] = sum(-1)^|S| E_a B_b.
Can we show that the Fourier transform hat{M}(j) = 0 for j >= 1
using the structure of E_a and B_b under cyclic shift?

For circulant T with shift sigma:
  E_{sigma(a)}(sigma(S)+sigma(a)) = E_a(S+a)
  B_{sigma(b)}(sigma(R)+sigma(b)) = B_b(R+b)

This gives M[sigma(a), sigma(b)] = M[a,b] (circulant property).
But does NOT force M[0,k] = 0.

We need ADDITIONAL structure. Candidates:
1. The inclusion-exclusion sign (-1)^|S| interacts with the cyclic structure
2. Path reversal (B_b = signed E_b) under converse
3. The tournament constraint T[i,j] + T[j,i] = 1

Let me compute the Fourier modes of E and B separately.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np

def make_circulant_tournament(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def E_v(T, verts, v):
    """Count Hamiltonian paths in T[verts] ending at v."""
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def B_v(T, verts, v):
    """Count Hamiltonian paths in T[verts] starting at v."""
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[0] != v:
            continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            val += sign * E_v(T, S_set, a) * B_v(T, R_set, b)
        return val

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count


# ============================================================
# Fourier analysis of E and B for circulant tournaments
# ============================================================
print("=" * 70)
print("Fourier analysis of E_a(S+a) and B_b(R+b)")
print("=" * 70)

n = 5
S_gen = {1, 2}  # Paley T_5
T = make_circulant_tournament(n, S_gen)
H = count_H(T, n)

print(f"\nPaley T_5 (S={{1,2}}): H={H}")

# For each subset size s, compute E_a(S+a) for a=0,...,n-1
# and take the DFT to see Fourier structure
for s_size in range(n):
    # All subsets of size s_size from [n]\{a}
    print(f"\n  Subset size |S|={s_size} (paths on {s_size+1} vertices):")

    # For each a, sum E_a over all subsets of size s_size
    E_sums = []
    B_sums = []
    for a in range(n):
        others = [v for v in range(n) if v != a]
        e_sum = 0
        b_sum = 0
        for S_tuple in combinations(others, s_size):
            S_set = sorted(list(S_tuple) + [a])
            e_sum += E_v(T, S_set, a)
            b_sum += B_v(T, S_set, a)
        E_sums.append(e_sum)
        B_sums.append(b_sum)

    print(f"    E_a sums: {E_sums}")
    print(f"    B_a sums: {B_sums}")

    # DFT
    E_dft = np.fft.fft(np.array(E_sums, dtype=complex))
    B_dft = np.fft.fft(np.array(B_sums, dtype=complex))
    print(f"    DFT(E): {[f'{d.real:.2f}+{d.imag:.2f}i' for d in E_dft]}")
    print(f"    DFT(B): {[f'{d.real:.2f}+{d.imag:.2f}i' for d in B_dft]}")


# ============================================================
# For the off-diagonal, we need the PAIRED Fourier structure
# ============================================================
print("\n" + "=" * 70)
print("Paired Fourier: E_a(S+a) * B_b(R+b) structure")
print("=" * 70)

# M[0,k] = sum_{S subset U_k} (-1)^|S| E_0(S+0) B_k(R+k)
# where U_k = [n]\{0,k}

# Under cyclic shift sigma: M[0,k] -> M[1,k+1]
# We want: for EACH k != 0, sum_S (-1)^|S| E_0(S+0) B_k(R+k) = 0

# Key observation: the U set changes with k.
# U_1 = {2,3,4}, U_2 = {1,3,4}, U_3 = {1,2,4}, U_4 = {1,2,3}

# Maybe we can use the converse identity:
# E_a(S+a, T) = B_a(S+a, T^op) * (-1)^|S|   (path reversal)
# T^op for Paley T_5: T^op has generating set {3,4} = {n-s for s in {1,2}}

T_op = make_circulant_tournament(n, {n-s for s in S_gen})
H_op = count_H(T_op, n)
print(f"\n  T_op generating set: {sorted({n-s for s in S_gen})}, H(T^op) = {H_op}")

# Verify path reversal: E_a(V, T) = B_a(V, T^op) for Hamiltonian paths
for a in range(n):
    ea_T = E_v(T, list(range(n)), a)
    ba_Top = B_v(T_op, list(range(n)), a)
    print(f"  E_{a}([n], T) = {ea_T}, B_{a}([n], T^op) = {ba_Top}", end="")
    # Actually path reversal gives: if P is path in T, rev(P) is path in T^op
    # P ends at a => rev(P) starts at a
    # So E_a(V, T) = B_a(V, T^op)
    print(f"  {'MATCH' if ea_T == ba_Top else 'MISMATCH'}")


# ============================================================
# Connection to position matrix and signed position
# ============================================================
print("\n" + "=" * 70)
print("Position matrix and transfer matrix connection")
print("=" * 70)

# P[v,k] = # Ham paths where v is at position k (0-indexed)
# M[a,a] = sum_P (-1)^{pos(a,P)} = sum_k (-1)^k P[a,k]
# This is the ALTERNATING sum of position counts.

# For M[a,b] (a != b), can we express it in terms of positions?
# M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
# E_a(S+a) counts paths in S+a ending at a
# B_b(R+b) counts paths in R+b starting at b

# A Ham path P decomposes as: P = path in S+a ending at a, then a->b edge(?), then path in R+b starting at b
# NO! There's no a->b edge in general. The M entry comes from inclusion-exclusion,
# not from actual path concatenation.

# But! Consider the "split" at vertex positions:
# If we split a Ham path P at position j (so P[j] = a):
#   prefix P[0:j+1] ends at a, has j+1 vertices
#   suffix P[j:n] starts at a, has n-j vertices
# This is NOT the same as E_a * B_b.

# However, if we condition on (position of a = j, position of b = k):
#   The # of such paths is related to E_a and B_b counts on the RIGHT subsets.

# Let's define P2[a,j,b,k] = # Ham paths where a is at position j and b is at position k.
P2 = np.zeros((n, n, n, n), dtype=int)
for perm in permutations(range(n)):
    prod = 1
    for k in range(n-1):
        prod *= T.get((perm[k], perm[k+1]), 0)
    if prod > 0:
        for i in range(n):
            for j in range(n):
                if i != j:
                    pos_i = list(perm).index(i)
                    pos_j = list(perm).index(j)
                    P2[i, pos_i, j, pos_j] += 1

# M[a,b] should be related to sum_{j,k} (-1)^? P2[a,j,b,k]
# Let's check: does sum_{j,k} (-1)^j P2[a,j,b,k] = ?
# Or: sum_{j,k} (-1)^{j+k} P2[a,j,b,k] = ?

for a in range(n):
    for b in range(a+1, n):
        M_val = compute_M_entry(T, n, a, b)

        # Various signed sums over P2
        s1 = sum((-1)**j * P2[a,j,b,k] for j in range(n) for k in range(n))
        s2 = sum((-1)**k * P2[a,j,b,k] for j in range(n) for k in range(n))
        s3 = sum((-1)**(j+k) * P2[a,j,b,k] for j in range(n) for k in range(n))
        s4 = sum((-1)**(j-k) * P2[a,j,b,k] for j in range(n) for k in range(n))
        s5 = sum((-1)**j * (-1)**(k>j) * P2[a,j,b,k] for j in range(n) for k in range(n))

        # What about: sum where j < k weighted by (-1)^j - sum where j > k weighted by (-1)^j?
        s6 = sum((-1)**j * P2[a,j,b,k] for j in range(n) for k in range(j+1, n))
        s7 = sum((-1)**j * P2[a,j,b,k] for j in range(n) for k in range(j))

        matches = []
        for name, val in [('(-1)^j', s1), ('(-1)^k', s2), ('(-1)^(j+k)', s3),
                          ('(-1)^(j-k)', s4), ('(-1)^j*(j<k)', s6),
                          ('(-1)^j*(j>k)', s7), ('s6-s7', s6-s7)]:
            if val == M_val:
                matches.append(name)

        if a <= 1 and b <= 2:
            print(f"  P2[{a},*,{b},*]: M={M_val}, (-1)^j={s1}, (-1)^k={s2}, "
                  f"(-1)^(j+k)={s3}, (-1)^j(j<k)={s6}, (-1)^j(j>k)={s7}")
            if matches:
                print(f"    MATCHES: {matches}")


# ============================================================
# ANOTHER APPROACH: Express M[a,b] directly via position matrix
# ============================================================
print("\n" + "=" * 70)
print("Direct formula: M[a,b] from position-pair matrix")
print("=" * 70)

# From the definition, M[a,a] = sum_P (-1)^{pos(a,P)}
# What about M[a,b]?
#
# Recall: M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
# S is a subset of U = [n]\{a,b}
# S+a is a subset containing a (and |S| other vertices)
# R+b is the complement containing b
#
# Each Ham path P of T can be "split" by which vertices are in S+a vs R+b.
# For a given P, vertex v is in S+a if v is among the first few vertices of P
# (those between position 0 and some cutpoint), or something else?
#
# Actually no, S is summed over ALL subsets, not determined by P.
#
# Let me try the ALTERNATING SUM interpretation differently.
# Define: f(a, V) = sum of Hamiltonian paths in T[V] where a is at EVEN position
#                 - sum where a is at ODD position
# Then M[a,a] = f(a, [n]).
#
# For M[a,b], the inclusion-exclusion formula sums over all ways to split
# [n]\{a,b} into S and R, with (-1)^|S| * (paths ending at a in S+a) * (paths starting at b in R+b).
#
# What if we interpret this as: enumerate pairs (path ending at a, path starting at b)
# on complementary vertex sets, with sign?

# Key insight: this is a CONVOLUTION. If we define
#   alpha(W, a) = E_a(W)  for a in W
#   beta(W, b)  = B_b(W)  for b in W
# then M[a,b] = sum_{S+a + R+b = [n]} (-1)^|S| alpha(S+a, a) beta(R+b, b)
#
# This is the "subset convolution" of alpha and beta.
# For circulant tournaments, alpha and beta have cyclic symmetry.

# Let me compute: for each vertex a, what is the total
# sum of (-1)^|S| * E_a(S+a) over ALL subsets S of [n]\{a}?
print(f"\n  Total signed E_a and B_a:")
for a in range(n):
    total_signed_E = 0
    total_signed_B = 0
    others = [v for v in range(n) if v != a]
    for mask in range(1 << len(others)):
        S_list = [others[k] for k in range(len(others)) if mask & (1 << k)]
        S_set = sorted(set(S_list) | {a})
        sz = len(S_list)
        ea = E_v(T, S_set, a)
        ba = B_v(T, S_set, a)
        total_signed_E += (-1)**sz * ea
        total_signed_B += (-1)**sz * ba
    print(f"    a={a}: sum (-1)^|S| E_a = {total_signed_E}, sum (-1)^|S| B_a = {total_signed_B}")


# ============================================================
# The alternating sum of E_a might have nice properties!
# ============================================================
print("\n" + "=" * 70)
print("Signed E_a sums by subset size")
print("=" * 70)

for a in range(n):
    others = [v for v in range(n) if v != a]
    by_size = {}
    for mask in range(1 << len(others)):
        S_list = [others[k] for k in range(len(others)) if mask & (1 << k)]
        S_set = sorted(set(S_list) | {a})
        sz = len(S_list)
        ea = E_v(T, S_set, a)
        if sz not in by_size:
            by_size[sz] = 0
        by_size[sz] += ea  # unsigned sum at each level
    layers = [by_size.get(k, 0) for k in range(n)]
    alt_layers = [(-1)**k * by_size.get(k, 0) for k in range(n)]
    print(f"  a={a}: E_a sums by |S|: {layers}, alternating: {alt_layers}, total alt: {sum(alt_layers)}")


# ============================================================
# Check if M[a,b] can be expressed as position-pair formula
# ============================================================
print("\n" + "=" * 70)
print("Testing position-pair formulas for M[a,b]")
print("=" * 70)

# Hypothesis: M[a,b] = sum_P (-1)^{pos(a,P)} * sign(pos(a,P), pos(b,P))
# where sign depends on relative positions of a and b.
#
# From the decomposition M[a,b] = sum_S (-1)^|S| E_a B_b, and
# the fact that S has |S| vertices besides a while R has |R| = n-2-|S| vertices besides b,
# and E_a is about paths ending at a in an (|S|+1)-vertex tournament while
# B_b is about paths starting at b in an (n-1-|S|)-vertex tournament,
# the "position" interpretation would be: pos(a,P) = |S| = position of a counting from the right.

# Let me check: M[a,a] = sum_P (-1)^{pos(a,P)}. Here pos is 0-indexed from left.
# Can we write M[a,b] = sum_P (-1)^{pos(a,P)} for paths where b has specific position?

# Actually, let me compute sum_P (-1)^{pos(a,P)} * [b at position k]
for a, b in [(0,1), (0,2)]:
    M_val = compute_M_entry(T, n, a, b)
    print(f"\n  M[{a},{b}] = {M_val}")

    for k in range(n):
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for q in range(n-1):
                prod *= T.get((perm[q], perm[q+1]), 0)
            if prod > 0:
                pos_a = list(perm).index(a)
                pos_b = list(perm).index(b)
                if pos_b == k:
                    val += (-1)**pos_a
        print(f"    sum_P (-1)^pos(a) [pos(b)={k}] = {val}")

    # What about sum_P (-1)^{pos(a,P)+pos(b,P)}?
    val_sum = 0
    for perm in permutations(range(n)):
        prod = 1
        for q in range(n-1):
            prod *= T.get((perm[q], perm[q+1]), 0)
        if prod > 0:
            pos_a = list(perm).index(a)
            pos_b = list(perm).index(b)
            val_sum += (-1)**(pos_a + pos_b)
    print(f"    sum_P (-1)^(pos(a)+pos(b)) = {val_sum}")

    # What about sum_P (-1)^{pos(a)} * (-1)^{pos(b)}?
    # Same as above

    # What about sum_P (-1)^{min(pos(a),pos(b))}?
    val_min = 0
    val_max = 0
    for perm in permutations(range(n)):
        prod = 1
        for q in range(n-1):
            prod *= T.get((perm[q], perm[q+1]), 0)
        if prod > 0:
            pos_a = list(perm).index(a)
            pos_b = list(perm).index(b)
            val_min += (-1)**min(pos_a, pos_b)
            val_max += (-1)**max(pos_a, pos_b)
    print(f"    sum_P (-1)^min(pos) = {val_min}, (-1)^max(pos) = {val_max}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
