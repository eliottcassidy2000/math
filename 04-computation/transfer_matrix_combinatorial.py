#!/usr/bin/env python3
"""
COMBINATORIAL MEANING of the transfer matrix M[a,b].

M[a,b] = sum_{S subset U} (-1)^|S| E_a(S+a) B_b(R+b)

where U = V\{a,b}, R = U\S.

QUESTION 1: What does M[a,a] count?
  M[a,a] = sum_perm (-1)^{pos_of_a} * [perm is Ham path]
  This is the SIGNED position count of vertex a in Hamiltonian paths.

QUESTION 2: What does M[a,b] (a != b) represent?
  It's an inclusion-exclusion over internal vertices.
  When all terms are positive (no cancellation), it counts paths from a to b.
  The alternating sign creates cancellations.

QUESTION 3: How does M relate to the Redei-Berge symmetric function U_T?
  The Schur expansion [s_{(k,1^{n-k})}]U_T counts paths with descent at position k.
  Is M[a,b] a refinement of this?

kind-pasteur-2026-03-06-S25 (continuation)
"""

from itertools import permutations
import numpy as np

def make_circulant(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count


# ============================================================
# QUESTION 1: M[a,a] = signed position sum
# ============================================================
print("=" * 70)
print("M[a,a]: signed position of vertex a in Hamiltonian paths")
print("=" * 70)

for n, name, S_set in [(3, "T_3 cyclic", {1}), (5, "Paley T_5", {1,2}),
                         (5, "circ {1,3}", {1,3})]:
    T = make_circulant(n, S_set)

    for a in range(n):
        # M[a,a] = sum over Ham paths of (-1)^{position of a}
        val = 0
        total_paths = 0
        pos_counts = [0]*n
        for perm in permutations(range(n)):
            prod = 1
            for k in range(len(perm)-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                total_paths += 1
                pos = list(perm).index(a)
                val += (-1)**pos
                pos_counts[pos] += 1

        print(f"\n  {name}, vertex {a}: M[{a},{a}] = {val}")
        print(f"    H(T) = {total_paths}")
        print(f"    Position distribution of {a}: {pos_counts}")
        print(f"    Signed sum: {' + '.join(f'({(-1)**k})*{pos_counts[k]}' for k in range(n))} = {val}")

        # Only print for vertex 0 to avoid clutter
        if a == 0:
            continue
        break


# ============================================================
# For VT tournaments: M[a,a] = H/n for all a
# ============================================================
print("\n" + "=" * 70)
print("Vertex-transitive: M[a,a] = H/n => position balance")
print("=" * 70)

print("""
If T is vertex-transitive, then every vertex has the same position distribution.
M[a,a] = H/n means: sum_k (-1)^k * P_k = H/n
where P_k = number of Ham paths with vertex a at position k.

Since sum_k P_k = H (each path counted once), we have:
  sum_k P_k = H
  sum_k (-1)^k P_k = H/n

This gives:
  sum_{even k} P_k = H(n+1)/(2n)
  sum_{odd k} P_k = H(n-1)/(2n)

So vertex-transitivity creates a specific even/odd position imbalance!
""")

# Verify
for n, name, S_set in [(3, "T_3", {1}), (5, "Paley T_5", {1,2}), (7, "Paley T_7", {1,2,4})]:
    T = make_circulant(n, S_set)
    H = count_H(T, n)

    # Compute position distribution for vertex 0
    pos_counts = [0]*n
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            pos = list(perm).index(0)
            pos_counts[pos] += 1

    even_sum = sum(pos_counts[k] for k in range(n) if k % 2 == 0)
    odd_sum = sum(pos_counts[k] for k in range(n) if k % 2 == 1)
    signed_sum = sum((-1)**k * pos_counts[k] for k in range(n))

    print(f"\n  {name} (n={n}): H = {H}, H/n = {H/n:.1f}")
    print(f"    Position dist of v=0: {pos_counts}")
    print(f"    Even positions: {even_sum}, Odd positions: {odd_sum}")
    print(f"    Signed sum: {signed_sum} (expected H/n = {H/n:.1f})")
    print(f"    Expected even: H(n+1)/(2n) = {H*(n+1)/(2*n):.1f}")
    print(f"    Expected odd:  H(n-1)/(2n) = {H*(n-1)/(2*n):.1f}")


# ============================================================
# QUESTION 2: M[a,b] for a != b - what does it mean?
# ============================================================
print("\n" + "=" * 70)
print("M[a,b] (a != b): inclusion-exclusion decomposition")
print("=" * 70)

# For n=3, M = I, so M[0,1] = 0. Let's see why.
n = 3
T = make_circulant(n, {1})
print(f"\n  n=3, cyclic 0->1->2->0:")
print(f"  Ham paths: 0->1->2, 2->0->1, 1->2->0")
print(f"  For M[0,1]: a=0, b=1, U={{2}}")
print(f"    S=empty: E_0({{0}}) * B_1({{1,2}}) = 1 * B_1({{1,2}})")
B12 = T.get((1,2), 0)  # path 1->2
print(f"      B_1({{1,2}}) = t(1,2) = {B12}")
print(f"    S={{2}}: (-1)^1 * E_0({{0,2}}) * B_1({{1}})")
E02 = T.get((2,0), 0)  # path 2->0 (ending at 0)
print(f"      E_0({{0,2}}) = t(2,0) = {E02}")
print(f"      B_1({{1}}) = 1")
print(f"    M[0,1] = {B12} + (-1)*{E02}*1 = {B12 - E02}")

print(f"\n  For M[0,2]: a=0, b=2, U={{1}}")
print(f"    S=empty: E_0({{0}}) * B_2({{1,2}}) = 1 * t(2,1)")
B21 = T.get((2,1), 0)
print(f"      B_2({{1,2}}) = t(2,1) = {B21}")
print(f"    S={{1}}: (-1)^1 * E_0({{0,1}}) * B_2({{2}})")
E01 = T.get((1,0), 0)
print(f"      E_0({{0,1}}) = t(1,0) = {E01}")
print(f"    M[0,2] = {B21} + (-1)*{E01}*1 = {B21 - E01}")


# ============================================================
# QUESTION 3: Deeper structure - M as alternating path decomposition
# ============================================================
print("\n" + "=" * 70)
print("DEEP STRUCTURE: M[a,b] as signed split-path count")
print("=" * 70)

print("""
M[a,b] = sum_S (-1)^|S| * E_a(S+a) * B_b(R+b)

This is an inclusion-exclusion that counts "split paths":
  - Left half: a Ham path through S+{a} ending at a
  - Right half: a Ham path through R+{b} starting at b
  - These two paths together cover all vertices but DON'T form a Ham path
    (there's no edge connecting them at the split point)
  - The inclusion-exclusion adjusts for overcounting

CLAIM: M[a,b] counts the NET contribution of vertex a's position and
vertex b's position to Hamiltonian path structure.

For a = b:
  M[a,a] = alternating position sum of a in Ham paths

For a != b:
  M[a,b] = 0 for VT tournaments at odd n (all info is on diagonal)
  M[a,b] measures "endpoint asymmetry" between a and b
""")

# Test: for NON-VT tournaments, what does M look like?
print("=" * 70)
print("NON-VT tournaments: transfer matrix structure")
print("=" * 70)

# Transitive tournament on 5 vertices
T = {}
for i in range(5):
    for j in range(5):
        if i != j:
            T[(i,j)] = 1 if i < j else 0

n = 5
H = count_H(T, n)
print(f"\n  Transitive T_5: H = {H}")

# Position distribution per vertex
for v in range(n):
    pos_counts = [0]*n
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            pos = list(perm).index(v)
            pos_counts[pos] += 1
    signed = sum((-1)**k * pos_counts[k] for k in range(n))
    print(f"    v={v}: positions={pos_counts}, signed_sum={signed}")

# M for transitive T_5
from itertools import permutations
def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(len(perm)-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod != 0:
                pos = list(perm).index(a)
                val += (-1)**pos * prod
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
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
                    prod *= T.get((p[k], p[k+1]), 0)
                ea += prod
            if len(S_set) == 1: ea = 1
            bb2 = 0
            for p in permutations(sorted(R_set)):
                if p[0] != b: continue
                prod = 1
                for k in range(len(p)-1):
                    prod *= T.get((p[k], p[k+1]), 0)
                bb2 += prod
            if len(R_set) == 1: bb2 = 1
            val += sign * ea * bb2
        return val

M = np.zeros((n,n), dtype=int)
for a in range(n):
    for b in range(n):
        M[a,b] = compute_M_entry(T, n, a, b)

print(f"\n    M (transitive T_5) =")
for row in M:
    print(f"      {list(row)}")
print(f"    trace = {int(np.trace(M))}")
print(f"    eigenvalues: {sorted(np.linalg.eigvalsh(M))[::-1]}")

# Also H(a->b) matrix
Hab = np.zeros((n,n), dtype=int)
for a in range(n):
    for b in range(n):
        if a == b: continue
        for perm in permutations(range(n)):
            if perm[0] != a or perm[-1] != b: continue
            prod = 1
            for k in range(len(perm)-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            Hab[a,b] += prod

print(f"\n    H(a->b) =")
for row in Hab:
    print(f"      {list(row)}")

# Key observation: for transitive tournament, only path is 0->1->2->3->4
# So Hab[0,4] = 1, all others 0
# M should reflect this single path


# ============================================================
# DISCOVERY: M eigenvalues and tournament structure
# ============================================================
print("\n" + "=" * 70)
print("M eigenvalues for different tournament types")
print("=" * 70)

import random
random.seed(42)

# Generate several n=5 tournaments of different types
tournaments = []

# 1. Cyclic (VT)
tournaments.append(("Cyclic {1,2}", make_circulant(5, {1,2})))
tournaments.append(("Cyclic {1,3}", make_circulant(5, {1,3})))

# 2. Transitive
T = {}
for i in range(5):
    for j in range(5):
        if i != j:
            T[(i,j)] = 1 if i < j else 0
tournaments.append(("Transitive", T))

# 3. Near-transitive (flip one arc)
T2 = dict(T)
T2[(0,4)] = 0; T2[(4,0)] = 1
tournaments.append(("Near-trans (flip 0->4)", T2))

# 4. Random regular
pairs = [(i,j) for i in range(5) for j in range(i+1, 5)]
for trial in range(1000):
    T3 = {}
    for (i,j) in pairs:
        if random.random() < 0.5:
            T3[(i,j)] = 1; T3[(j,i)] = 0
        else:
            T3[(i,j)] = 0; T3[(j,i)] = 1
    # Check regular
    if all(sum(T3.get((i,j),0) for j in range(5) if j != i) == 2 for i in range(5)):
        key = tuple(T3.get((i,j), 0) for (i,j) in pairs)
        is_cyclic = any(np.array_equal(
            np.array([[T3.get((i,j),0) for j in range(5)] for i in range(5)]),
            np.array([[make_circulant(5,S).get((i,j),0) for j in range(5)] for i in range(5)])
        ) for S in [{1,2},{1,3},{2,4},{3,4}])
        if not is_cyclic:
            tournaments.append(("Random regular", T3))
            break

n = 5
for name, T in tournaments:
    H = count_H(T, n)
    M = np.zeros((n,n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)
    evals = sorted(np.linalg.eigvalsh(M))[::-1]
    is_scalar = np.allclose(M, M[0,0]*np.eye(n))
    is_sym = np.allclose(M, M.T)
    print(f"\n  {name}: H={H}, symmetric={is_sym}")
    print(f"    M diagonal: {[int(M[i,i]) for i in range(n)]}")
    print(f"    Scalar: {is_scalar}")
    print(f"    Eigenvalues: {[round(e,2) for e in evals]}")
    if not is_scalar:
        print(f"    M =")
        for row in M:
            print(f"      {list(row)}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
