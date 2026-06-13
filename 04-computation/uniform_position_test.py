#!/usr/bin/env python3
"""
DISCOVERY: Paley tournaments have UNIFORM position distribution.

For Paley T_p: each vertex v appears at position k in exactly H/n Hamiltonian
paths, for every k = 0, 1, ..., n-1.

This is MUCH stronger than M[a,a] = H/n (which only says the signed sum
of position counts equals H/n).

QUESTION: Is uniform position distribution equivalent to M = (H/n)*I?
QUESTION: Do all VT tournaments have uniform position?
QUESTION: Do any non-VT regular tournaments have it?

kind-pasteur-2026-03-06-S25 (continuation)
"""

from itertools import permutations
import numpy as np
import random

def make_circulant(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def compute_position_matrix(T, n):
    """P[v,k] = number of Ham paths with vertex v at position k."""
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def is_regular(T, n):
    for i in range(n):
        if sum(T.get((i,j),0) for j in range(n) if j != i) != (n-1)//2:
            return False
    return True


# ============================================================
# Test VT tournaments at n=3, 5, 7
# ============================================================
print("=" * 70)
print("Position distribution for vertex-transitive tournaments")
print("=" * 70)

for n, name, S_set in [(3, "T_3", {1}),
                         (5, "Paley T_5 {1,2}", {1,2}),
                         (5, "circ {1,3}", {1,3}),
                         (7, "Paley T_7 {1,2,4}", {1,2,4}),
                         (7, "circ {1,3,5}", {1,3,5}),
                         (7, "circ {2,3,4}", {2,3,4})]:
    T = make_circulant(n, S_set)
    H = count_H(T, n)
    P = compute_position_matrix(T, n)

    # Check if position distribution is uniform
    is_uniform = np.all(P == H // n)

    print(f"\n  {name} (n={n}): H = {H}, H/n = {H/n:.1f}")
    if is_uniform:
        print(f"    Position distribution: UNIFORM ({H//n} at every position)")
    else:
        print(f"    Position distribution of v=0: {list(P[0])}")
        # Check if at least row-constant (same for each vertex)
        row_constant = all(np.array_equal(P[i], P[0]) for i in range(n))
        print(f"    Row-constant (same for all vertices): {row_constant}")
        if not row_constant:
            for v in range(min(n,4)):
                print(f"      v={v}: {list(P[v])}")


# ============================================================
# Test: ALL regular tournaments at n=5
# ============================================================
print("\n" + "=" * 70)
print("ALL regular n=5 tournaments: position distribution")
print("=" * 70)

random.seed(42)
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
seen = set()
results = []

for trial in range(100000):
    T = {}
    for (i,j) in pairs:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    if not is_regular(T, n):
        continue

    key = tuple(T.get((i,j),0) for (i,j) in pairs)
    if key in seen:
        continue
    seen.add(key)

    H = count_H(T, n)
    P = compute_position_matrix(T, n)
    is_uniform = np.all(P == H // n)

    # Check M = scalar
    # Actually since ALL n=5 regular have M=3I, they all have signed position = 3
    # But do they all have UNIFORM position?

    results.append((H, is_uniform, P))

print(f"  Found {len(results)} distinct regular n=5 tournaments")

uniform_count = sum(1 for _, u, _ in results if u)
print(f"  Uniform position: {uniform_count}/{len(results)}")

# Show non-uniform examples if any
for H, is_uniform, P in results:
    if not is_uniform:
        print(f"\n  NON-UNIFORM example: H={H}")
        for v in range(n):
            print(f"    v={v}: {list(P[v])}")
        break
else:
    print("  ALL regular n=5 tournaments have uniform position distribution!")


# ============================================================
# Test: Some regular tournaments at n=7
# ============================================================
print("\n" + "=" * 70)
print("Regular n=7 tournaments: position distribution")
print("=" * 70)

random.seed(123)
n = 7
pairs7 = [(i,j) for i in range(n) for j in range(i+1, n)]
seen7 = set()
results7 = []

for trial in range(500000):
    T = {}
    for (i,j) in pairs7:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    if not is_regular(T, n):
        continue

    key = tuple(T.get((i,j),0) for (i,j) in pairs7)
    if key in seen7:
        continue
    seen7.add(key)

    H = count_H(T, n)
    P = compute_position_matrix(T, n)
    is_uniform = np.all(P == H // n)

    results7.append((H, is_uniform, P))

    if len(results7) >= 20:
        break

print(f"  Found {len(results7)} distinct regular n=7 tournaments")

for H, is_uniform, P in results7:
    if is_uniform:
        print(f"  H={H}: UNIFORM ({H//n} at every position)")
    else:
        print(f"  H={H}: NON-UNIFORM")
        print(f"    v=0: {list(P[0])}")
        print(f"    v=1: {list(P[1])}")

uniform_count7 = sum(1 for _, u, _ in results7 if u)
print(f"\n  Uniform: {uniform_count7}/{len(results7)}")


# ============================================================
# KEY QUESTION: Is uniform position equivalent to M = scalar?
# ============================================================
print("\n" + "=" * 70)
print("Is uniform position <=> M = (H/n)*I?")
print("=" * 70)

print("""
Uniform position: P[v,k] = H/n for all v, k
This means: every vertex appears at every position equally often.

M[a,a] = sum_k (-1)^k P[a,k]
If P[a,k] = H/n for all k:
  M[a,a] = (H/n) * sum_k (-1)^k = (H/n) * [n is odd ? 1 : 0]
  For odd n: M[a,a] = H/n.

For odd n, M[a,a] = H/n is NECESSARY for M = (H/n)*I.
Uniform position IMPLIES M[a,a] = H/n.
But M[a,a] = H/n does NOT imply uniform position!

The signed position sum H/n is weaker than uniform position.
So uniform position is STRONGER than M = (H/n)*I.

QUESTION: Do we actually have uniform position for all regular n=5?
Answer: YES (verified above).

QUESTION: Is this a consequence of regularity at n=5?
At n=5, all regular tournaments are vertex-transitive (only 2 circulant
isomorphism classes). So uniform position follows from VT.

At n=7, there are non-VT regular tournaments. Do THOSE have uniform position?
""")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
