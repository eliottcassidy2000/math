#!/usr/bin/env python3
"""
Connection between u_T(m) (THM-077) and P(u,x) (THM-074).

THM-074: P(u,x) = P_n(u,0) + sum_I c_I * x^|I| * P_{n-S_I}(u,0) * (u-2)^{S_I/2}
THM-077: u_T(m) = sum_j sw(j) * m^{n-2j}

KEY OBSERVATION:
- In THM-074, the grouping is by S_I = sum(l_i - 1) = total "consumed" positions
- In THM-077, the grouping is by j = sum (|C_i|-1)/2 = S_I/2
- So sw(j) groups independent sets by S = 2j, same as THM-074!

The connection: u_T(m) is the evaluation of P(u,x) at... what?

Actually: G_T(t,x) = t^m * P(u,x) where u = t + 1/t.
And u_T(m) = ps_1(U_T)(m).

The key: ps_1(U_T)(m) evaluates U_T at p_k -> m for all k.
This is related to the Grinberg-Stanley symmetric function, not directly to P(u,x).

But the STRUCTURE is the same: both group by S_I.

Let me verify: does sw(j) equal the sum of c_I * 2^|I| over all I with S_I = 2j?

For a single cycle C of size l: S_C = l-1, |{C}| = 1, c_C = 1 (if C exists in Omega)
Weight = 2^1 = 2. And j = (l-1)/2.

For a pair {C_1, C_2} disjoint: S = (l_1-1)+(l_2-1), |I|=2, weight = 2^2 = 4.
j = S/2.

So sw(j) = sum_{I indep, S_I=2j} 2^{|I|}, and THM-074's grouping at u=0, x=2
would give sum_{I, S_I=2j} c_I * 2^|I| * P_{n-2j}(0,0) * (-2)^j.

At u=0 (i.e., t=i): P_k(0,0) = 2^{(k-1)/2} * T_k (tangent number).

So the total is:
P(0, 2) = sum_j [sum_{I, S_I=2j} c_I * 2^|I|] * P_{n-2j}(0,0) * (-2)^j
         = sum_j sw(j) * 2^{(n-2j-1)/2} * T_{n-2j} * (-2)^j

And P(0, 2) should equal... G_T(i, 2) / i^m.
G_T(t, 2) at t=i: this is the Hamiltonian path count at t=i...

Actually, from THM-074: G_T(0, x) = I(Omega(T), x), so G_T(0, 2) = H(T).
And P(u, 2) at u=0: P(0, 2) = G_T(i, 2) / i^m... but G_T(t,x) = t^m * P(t+1/t, x).

Hmm, actually G_T(0, 2) = 0^m * P(...) which is 0 for m >= 1. That's wrong.

Let me re-read: "G_T(t,x) = t^m * P(u,x)" with "G_T(0, x) = I(Omega(T), x)".
This means t^m * P(t + 1/t, x) at t=0. As t->0, u = t + 1/t -> infinity.
P(u, x) must have the right behavior as u->infinity.

Actually, G_T(t, x) is a polynomial in t (not a Laurent series). The relation
G_T(t, x) = t^m * P(t + 1/t, x) means P is a Laurent polynomial in t.
But P(u, x) is a polynomial in u, and u = t + 1/t is a Laurent polynomial.

So G_T(t, x) = t^m * P(t + 1/t, x) is the conversion.
At t = 0: G_T(0, x) = I(Omega, x). But t^m * P(t+1/t, x) -> 0 as t->0 (since
t^m -> 0 and P(t+1/t, x) grows like 1/t^{deg_u P}... the limit needs care).

Actually, the correct interpretation: G_T(t, x) is a polynomial in t of degree n-1.
Write G_T(t, x) = sum_{k=0}^{n-1} c_k(x) * t^k.
Then G_T(0, x) = c_0(x) = I(Omega, x). ← This is the CONSTANT TERM of G_T in t.

So the u_T(m) polynomial is u_T(m) = ps_1(U_T)(m), which is about the symmetric function,
not directly about G_T(t,x) or P(u,x).

But the GROUPING by S = 2j is the same in both! This is the key structural connection.

Let me just verify the sw decomposition numerically.

opus-2026-03-07-S39
"""
from itertools import permutations, combinations
from collections import defaultdict


def make_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for k, (i,j) in enumerate(edges):
        if bits & (1 << k):
            A[j][i] = 1
        else:
            A[i][j] = 1
    return A


def find_all_cycles(n, A):
    edge_set = {(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]}
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            seen = set()
            for p in permutations(verts):
                if all((p[i], p[(i+1)%length]) in edge_set for i in range(length)):
                    min_idx = list(p).index(min(p))
                    canon = tuple(list(p)[min_idx:] + list(p)[:min_idx])
                    if canon not in seen:
                        seen.add(canon)
                        cycles.append((frozenset(verts), length))
    return cycles


def compute_grouped_sums(n, A):
    """Group independent sets by (|I|, S_I) and by j = S_I/2."""
    cycles = find_all_cycles(n, A)
    nc = len(cycles)

    # Group by |I| (standard alpha) and by j = S/2 (size-weighted)
    by_size = defaultdict(int)  # alpha_k
    by_j = defaultdict(int)     # sw(j) = sum 2^|I| for S_I = 2j
    by_size_and_j = defaultdict(int)  # (k, j) -> count of I with |I|=k, S_I=2j

    def backtrack(idx, used_verts, k, j_sum):
        by_size[k] += 1
        by_j[j_sum] += 2**k
        by_size_and_j[(k, j_sum)] += 1
        for i in range(idx, nc):
            vset, length = cycles[i]
            if not (vset & used_verts):
                backtrack(i+1, used_verts | vset, k+1, j_sum + (length-1)//2)

    backtrack(0, frozenset(), 0, 0)
    return by_size, by_j, by_size_and_j


# === n=5 C_5 example ===
print("=== C_5 (Paley T_5) ===")
A = [[0]*5 for _ in range(5)]
for i in range(5):
    A[i][(i+1)%5] = 1
    A[i][(i+2)%5] = 1

alpha, sw, size_j = compute_grouped_sums(5, A)
cycles = find_all_cycles(5, A)
size_counts = defaultdict(int)
for _, l in cycles:
    size_counts[l] += 1

print(f"Cycles by size: {dict(sorted(size_counts.items()))}")
print(f"alpha (by |I|): {dict(sorted(alpha.items()))}")
print(f"sw (by j=S/2): {dict(sorted(sw.items()))}")
print(f"(|I|, j) counts: {dict(sorted(size_j.items()))}")
print()
print("THM-074 grouping (by S=2j):")
for j in sorted(sw.keys()):
    S = 2*j
    contrib_types = [(k, jj) for (k, jj), cnt in size_j.items() if jj == j and k > 0]
    print(f"  S={S} (j={j}): sw={sw[j]}, types={(sorted(contrib_types))}")

# === n=7 Paley T_7 ===
print("\n=== Paley T_7 ===")
A = [[0]*7 for _ in range(7)]
QR = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in QR:
            A[i][j] = 1

alpha, sw, size_j = compute_grouped_sums(7, A)
cycles = find_all_cycles(7, A)
size_counts = defaultdict(int)
for _, l in cycles:
    size_counts[l] += 1

print(f"Cycles by size: {dict(sorted(size_counts.items()))}")
print(f"alpha (by |I|): {dict(sorted(alpha.items()))}")
print(f"sw (by j=S/2): {dict(sorted(sw.items()))}")
print(f"H = I(Omega,2) = {sum(2**k * alpha[k] for k in alpha)}")
print(f"H from sw = {sum(sw.values())}")

# Show (|I|, j) decomposition
print(f"\n(|I|, j) decomposition:")
for (k, j), cnt in sorted(size_j.items()):
    if cnt > 0:
        S = 2*j
        cycle_info = f"S={S}"
        print(f"  |I|={k}, j={j} ({cycle_info}): {cnt} independent sets")

# What does this look like in u_T(m)?
print(f"\nu_T(m) = {' + '.join(f'{sw[j]}*m^{7-2*j}' for j in sorted(sw.keys()))}")
# Verify u_T(1) = H
print(f"u_T(1) = {sum(sw.values())} = H(T)")

# === Compare with some non-regular n=7 tournaments ===
print("\n=== Non-regular n=7 example ===")
import random
random.seed(42)
A = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(i+1, 7):
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

alpha, sw, size_j = compute_grouped_sums(7, A)
cycles = find_all_cycles(7, A)
size_counts = defaultdict(int)
for _, l in cycles:
    size_counts[l] += 1

print(f"Cycles by size: {dict(sorted(size_counts.items()))}")
print(f"alpha: {dict(sorted(alpha.items()))}")
print(f"sw: {dict(sorted(sw.items()))}")
print(f"H = {sum(2**k * alpha[k] for k in alpha)}")

print(f"\n(|I|, j) decomposition:")
for (k, j), cnt in sorted(size_j.items()):
    if cnt > 0 and k > 0:
        print(f"  |I|={k}, j={j} (S={2*j}): {cnt} independent sets, contrib to sw: {cnt * 2**k}")
