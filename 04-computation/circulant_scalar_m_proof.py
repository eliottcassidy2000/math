#!/usr/bin/env python3
"""
OBSERVATION: ALL circulant tournaments on 7 vertices have M = (H/n)*I.
Is this true for all circulant tournaments?

A circulant tournament on n vertices (n odd) with generating set S:
  i -> j iff (j-i) mod n in S, where S and {n-s: s in S} partition {1,...,n-1}.

The cyclic symmetry sigma(i) = i+1 mod n is an automorphism.
So M[sigma(a), sigma(b)] = M[a,b] for all a,b.
This means M is circulant! A circulant matrix that is also symmetric (THM-030)
and has M[0,0] = H/n (since sum of row = H and circulant rows sum same).

Wait: for a symmetric circulant matrix, the entries satisfy M[i,j] = f(|i-j| mod n).
And f(0) = M[0,0] = H/n (from trace = H).
If M = (H/n)*I, then f(k) = 0 for k != 0.

Is f(k) = 0 guaranteed for circulant tournaments?

By the Key Identity (THM-030), for VT at odd n: B_b = E_b for all b.
For circulant tournaments, all vertices are equivalent, so VT holds.
So B_b = H/n = E_b for all b.

But M[a,b] being zero for a != b requires more:
  M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b) = 0

The cyclic symmetry forces M circulant. So M[0,k] = M[0,1] for k=1,...,n-1?
No, that's not right either. Circulant means M[i,j] = g((j-i) mod n).

Let's check: is M[0,k] the same for all k != 0?

kind-pasteur-2026-03-06-S25b (continuation)
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
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    ea += prod
            bb2 = 0
            if len(R_set) == 1:
                bb2 = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    bb2 += prod
            val += sign * ea * bb2
        return val

def compute_M(T, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        M[a, a] = compute_M_entry(T, n, a, a)
        for b in range(a+1, n):
            M[a, b] = compute_M_entry(T, n, a, b)
            M[b, a] = M[a, b]
    return M

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count


# ============================================================
# Test: M circulant structure for circulant tournaments
# ============================================================
print("=" * 70)
print("Circulant tournaments: M structure")
print("=" * 70)

for n in [3, 5, 7]:
    print(f"\n  n={n}:")
    # Find all valid circulant generating sets
    gen_sets = []
    for S_tuple in combinations(range(1, n), (n-1)//2):
        S = set(S_tuple)
        complement = {(n - s) % n for s in S}
        if not (S & complement) and (S | complement) == set(range(1, n)):
            gen_sets.append(S)

    for S in gen_sets:
        T = make_circulant_tournament(n, S)
        H = count_H(T, n)
        M = compute_M(T, n)

        # Check if M is circulant
        is_circulant = True
        for i in range(n):
            for j in range(n):
                if M[i,j] != M[(i+1)%n, (j+1)%n]:
                    is_circulant = False
                    break

        # Check if M is scalar (H/n * I)
        scalar = H // n
        is_scalar = np.array_equal(M, scalar * np.eye(n, dtype=int))

        # Show M[0,:] to see the circulant pattern
        row0 = [M[0, j] for j in range(n)]
        print(f"    S={sorted(S)}: H={H}, H/n={H/n:.1f}, M[0,:]={row0}, circulant={is_circulant}, scalar={is_scalar}")


# ============================================================
# WHY is M scalar for circulant tournaments?
# ============================================================
print("\n" + "=" * 70)
print("WHY scalar? Analyzing M[0,k] for k != 0")
print("=" * 70)

# For circulant tournaments, M[0,k] = M[0,1] for all k (by automorphism group).
# Wait, that's only true if the aut group is transitive on PAIRS, not just vertices.
# The aut group Z/nZ sends (0,k) to (1,k+1) to (2,k+2), etc.
# So M[0,k] = M[1,k+1] = M[2,k+2] etc. (circulant).
# But M[0,1] vs M[0,2] could be different!

# Unless the automorphism group is 2-transitive...

n = 5
print(f"\n  n={n}:")
for S in [{1,2}, {1,3}]:
    T = make_circulant_tournament(n, S)
    H = count_H(T, n)
    M = compute_M(T, n)
    print(f"\n    S={sorted(S)}: H={H}")
    print(f"    M =")
    for row in M:
        print(f"      {list(row)}")

    # Check: M[0,k] for different k
    print(f"    M[0,k]: {[M[0,k] for k in range(n)]}")

    # The circulant structure: M[0,k] depends only on k
    # But for M to be scalar, ALL M[0,k] for k!=0 must be 0.
    # This means: M[0,1] = 0, M[0,2] = 0, etc.

n = 7
print(f"\n  n={n}:")
for S in [{1,2,4}, {1,2,3}, {1,3,5}]:
    T = make_circulant_tournament(n, S)
    H = count_H(T, n)
    M = compute_M(T, n)
    print(f"\n    S={sorted(S)}: H={H}")
    print(f"    M[0,k]: {[M[0,k] for k in range(n)]}")


# ============================================================
# Non-circulant regular tournaments: is M still scalar?
# ============================================================
print("\n" + "=" * 70)
print("Non-circulant regular tournaments at n=7")
print("=" * 70)

# At n=7, are there regular (3-regular) tournaments that are NOT circulant?
import random
random.seed(42)

pairs = [(i,j) for i in range(7) for j in range(i+1, 7)]
found = 0

for trial in range(100000):
    T = {}
    for (i,j) in pairs:
        if random.random() < 0.5:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1

    # Check regular
    scores = tuple(sum(T.get((i,j),0) for j in range(7) if j != i) for i in range(7))
    if not all(s == 3 for s in scores):
        continue

    # Check if it's circulant
    is_circ = False
    for S_tuple in combinations(range(1, 7), 3):
        S = set(S_tuple)
        match = True
        for i in range(7):
            for j in range(7):
                if i != j:
                    expected = 1 if ((j-i) % 7) in S else 0
                    if T.get((i,j), 0) != expected:
                        match = False
                        break
            if not match:
                break
        if match:
            is_circ = True
            break

    if not is_circ:
        H = count_H(T, 7)
        M = compute_M(T, 7)
        scalar = H // 7
        is_scalar = np.array_equal(M, scalar * np.eye(7, dtype=int))

        print(f"  Non-circulant regular: H={H}, M scalar={is_scalar}")
        if not is_scalar:
            print(f"    M[0,:] = {[M[0,k] for k in range(7)]}")

        found += 1
        if found >= 5:
            break

if found == 0:
    print("  No non-circulant regular tournaments found (all regular n=7 are circulant?)")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
