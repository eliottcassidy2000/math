#!/usr/bin/env python3
"""
What exactly characterizes scalar M at odd n?

Known:
- VT => scalar M (THM-052, proved)
- scalar M <=> position-uniform (at n=3,5)
- scalar M =/=> VT (disproved at n=5)
- All scalar-M at n=5 are self-complementary (kind-pasteur verified)

The non-VT scalar-M tournament at n=5 has:
  Vertex orbits under Aut: {0}, {1,2,3}, {4}
  Anti-automorphisms swap 0 <-> 4
  => Aut(T) union Anti(T) acts TRANSITIVELY on V!

HYPOTHESIS: T has scalar M iff Aut(T) union Anti(T) acts transitively on V.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np

def transfer_matrix(A):
    n = len(A)
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_subset(A, S_verts, end=a)
                    bb = count_paths_subset(A, R_verts, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        if all(A[p[i]][p[i+1]] == 1 for i in range(len(p)-1)):
            count += 1
    return count

def automorphisms(A):
    n = len(A)
    return [perm for perm in permutations(range(n))
            if all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n))]

def anti_automorphisms(A):
    """Anti-automorphism: sigma such that A[sigma(i)][sigma(j)] = A[j][i] for all i,j.
    Equivalently: sigma is an isomorphism from T to T^op."""
    n = len(A)
    return [perm for perm in permutations(range(n))
            if all(A[perm[i]][perm[j]] == A[j][i] for i in range(n) for j in range(n))]

def orbits_of_perms(perms, n):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
    for perm in perms:
        for i in range(n):
            union(i, perm[i])
    orbs = {}
    for i in range(n):
        r = find(i)
        if r not in orbs:
            orbs[r] = []
        orbs[r].append(i)
    return list(orbs.values())

# =====================================================================
print("=" * 70)
print("n=5: Aut(T) UNION Anti(T) TRANSITIVE <=> SCALAR M?")
print("=" * 70)

n = 5
counts = {"scalar_transitive": 0, "scalar_not_transitive": 0,
          "nonscalar_transitive": 0, "nonscalar_not_transitive": 0}

for bits in range(1024):
    A = [[0]*n for _ in range(n)]
    edge_idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> edge_idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            edge_idx += 1

    H = sum(1 for p in permutations(range(n))
            if all(A[p[k]][p[k+1]] == 1 for k in range(n-1)))

    # Scalar M check
    if H % n != 0:
        scalar = False
    else:
        M = transfer_matrix(A)
        scalar = all(M[i][j] == (H//n if i == j else 0) for i in range(n) for j in range(n))

    # Aut + Anti transitivity check
    auts = automorphisms(A)
    antis = anti_automorphisms(A)
    all_perms = auts + antis
    orbs = orbits_of_perms(all_perms, n)
    transitive = len(orbs) == 1

    if scalar and transitive:
        counts["scalar_transitive"] += 1
    elif scalar and not transitive:
        counts["scalar_not_transitive"] += 1
        scores = tuple(sorted(sum(row) for row in A))
        print(f"  COUNTEREX (scalar, not trans): bits={bits}, H={H}, scores={scores}, "
              f"|Aut|={len(auts)}, |Anti|={len(antis)}, orbits={orbs}")
    elif not scalar and transitive:
        counts["nonscalar_transitive"] += 1
        scores = tuple(sorted(sum(row) for row in A))
        if counts["nonscalar_transitive"] <= 5:
            print(f"  COUNTEREX (non-scalar, trans): bits={bits}, H={H}, scores={scores}, "
                  f"|Aut|={len(auts)}, |Anti|={len(antis)}")
    else:
        counts["nonscalar_not_transitive"] += 1

print(f"\n  Results:")
print(f"    Scalar M + Aut∪Anti transitive: {counts['scalar_transitive']}")
print(f"    Scalar M + NOT transitive: {counts['scalar_not_transitive']}")
print(f"    Non-scalar + Aut∪Anti transitive: {counts['nonscalar_transitive']}")
print(f"    Non-scalar + NOT transitive: {counts['nonscalar_not_transitive']}")

if counts['scalar_not_transitive'] == 0 and counts['nonscalar_transitive'] == 0:
    print("\n  THEOREM: At n=5, scalar M <=> Aut(T) ∪ Anti(T) acts transitively on V")
elif counts['scalar_not_transitive'] == 0:
    print(f"\n  Scalar M => Aut∪Anti transitive (but {counts['nonscalar_transitive']} non-scalar are also transitive)")
elif counts['nonscalar_transitive'] == 0:
    print(f"\n  Aut∪Anti transitive => scalar M (but {counts['scalar_not_transitive']} scalar ones are not transitive)")

# =====================================================================
# n=3 check
# =====================================================================
print("\n" + "=" * 70)
print("n=3: SAME CHECK")
print("=" * 70)

n = 3
counts3 = {"scalar_transitive": 0, "scalar_not_transitive": 0,
           "nonscalar_transitive": 0, "nonscalar_not_transitive": 0}

for bits in range(8):
    A = [[0]*3 for _ in range(3)]
    pairs = [(0,1), (0,2), (1,2)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    H = sum(1 for p in permutations(range(3)) if all(A[p[k]][p[k+1]] == 1 for k in range(2)))
    M = transfer_matrix(A)
    scalar = all(M[i][j] == (H//3 if i == j else 0) for i in range(3) for j in range(3))

    auts = automorphisms(A)
    antis = anti_automorphisms(A)
    all_perms = auts + antis
    orbs = orbits_of_perms(all_perms, 3)
    transitive = len(orbs) == 1

    if scalar and transitive:
        counts3["scalar_transitive"] += 1
    elif scalar and not transitive:
        counts3["scalar_not_transitive"] += 1
    elif not scalar and transitive:
        counts3["nonscalar_transitive"] += 1
    else:
        counts3["nonscalar_not_transitive"] += 1

print(f"  Results:")
for k, v in counts3.items():
    print(f"    {k}: {v}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
