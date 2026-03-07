#!/usr/bin/env python3
"""
Verify the "NOT self-converse" finding for F_21 Cayley tournament.

The backtracking search said NO anti-automorphism exists for a specific
non-normal connection set. But the hash invariant matched T and T^op.

Let's verify by:
1. Double-checking the backtracking search correctness
2. Using a different isomorphism approach
3. Testing with a KNOWN self-converse case first

kind-pasteur-2026-03-06-S25e
"""

POW2 = [1, 2, 4]

def mul(a, b):
    return ((a[0] + POW2[a[1]] * b[0]) % 7, (a[1] + b[1]) % 3)

def inv(a):
    j_inv = (-a[1]) % 3
    i_inv = (-POW2[j_inv] * a[0]) % 7
    return (i_inv, j_inv)

elements = [(i, j) for j in range(3) for i in range(7)]
idx = {e: k for k, e in enumerate(elements)}
n = 21

nonid = [e for e in elements if e != (0, 0)]
pairs = []
seen = set()
for g in nonid:
    if g not in seen:
        pairs.append((g, inv(g)))
        seen.add(g)
        seen.add(inv(g))

def adj_matrix(S_set):
    A = [[0]*n for _ in range(n)]
    for a in elements:
        for b in elements:
            if a == b:
                continue
            diff = mul(inv(a), b)
            if diff in S_set:
                A[idx[a]][idx[b]] = 1
    return A

# ============================================================
# TEST 1: Verify backtracking on a KNOWN self-converse tournament
# ============================================================
print("=" * 70)
print("TEST 1: Known self-converse (normal S)")
print("=" * 70)

# Normal S = C1 union C3 = {(1,0),(2,0),(4,0)} union {(a,1): a in Z_7}
C1 = {(1,0), (2,0), (4,0)}
C3 = {(a, 1) for a in range(7)}
S_normal = C1 | C3
print(f"  S = {sorted(S_normal)}, |S| = {len(S_normal)}")

A_normal = adj_matrix(S_normal)

# Find anti-automorphism by backtracking
def backtrack_anti_auto(A, assignment, remaining, n):
    if not remaining:
        return assignment.copy()

    v = remaining[0]
    rest = remaining[1:]
    used = set(assignment.values())

    for target in range(n):
        if target in used:
            continue
        ok = True
        for u, su in assignment.items():
            if A[su][target] != A[v][u]:
                ok = False
                break
            if A[target][su] != A[u][v]:
                ok = False
                break
        if ok:
            assignment[v] = target
            result = backtrack_anti_auto(A, assignment, rest, n)
            if result is not None:
                return result
            del assignment[v]

    return None

# Order vertices by degree of constraint (most constrained first)
# For VT, all vertices are equivalent, so order doesn't matter much.
# But let's try vertex 0 first, then its out-neighbors, etc.
ordering = list(range(n))

result = backtrack_anti_auto(A_normal, {}, ordering, n)
if result is not None:
    sigma = [result[i] for i in range(n)]
    ok = all(A_normal[sigma[i]][sigma[j]] == A_normal[j][i]
             for i in range(n) for j in range(n) if i != j)
    print(f"  Anti-automorphism found! Verified: {ok}")
else:
    print(f"  ERROR: No anti-automorphism found for KNOWN self-converse tournament!")
    print(f"  This indicates a BUG in the backtracking search!")

# ============================================================
# TEST 2: The non-normal S that reportedly failed
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: Non-normal S (claimed NOT self-converse)")
print("=" * 70)

S_nonnormal = set()
for k, (g, gi) in enumerate(pairs):
    if k < 5:
        S_nonnormal.add(g)
    else:
        S_nonnormal.add(gi)

print(f"  S = {sorted(S_nonnormal)}, |S| = {len(S_nonnormal)}")

A_nonnormal = adj_matrix(S_nonnormal)

# Verify it's a tournament
for i in range(n):
    for j in range(i+1, n):
        assert A_nonnormal[i][j] + A_nonnormal[j][i] == 1

# Compare some structural invariants of T and T^op
A_op = [[A_nonnormal[j][i] for j in range(n)] for i in range(n)]

# Count 3-cycles
tc_T = sum(1 for i in range(n) for j in range(n) for k in range(n)
           if A_nonnormal[i][j] == 1 and A_nonnormal[j][k] == 1 and A_nonnormal[k][i] == 1)
tc_Top = sum(1 for i in range(n) for j in range(n) for k in range(n)
            if A_op[i][j] == 1 and A_op[j][k] == 1 and A_op[k][i] == 1)
print(f"  3-cycles in T: {tc_T}, in T^op: {tc_Top}")

# Count 4-cycles
fc_T = sum(1 for i in range(n) for j in range(n) for k in range(n) for l in range(n)
           if A_nonnormal[i][j]==1 and A_nonnormal[j][k]==1 and A_nonnormal[k][l]==1 and A_nonnormal[l][i]==1)
fc_Top = sum(1 for i in range(n) for j in range(n) for k in range(n) for l in range(n)
            if A_op[i][j]==1 and A_op[j][k]==1 and A_op[k][l]==1 and A_op[l][i]==1)
print(f"  4-cycles in T: {fc_T}, in T^op: {fc_Top}")

# For per-vertex: out-neighborhood edge count profile
def out_edge_profile(A):
    profiles = []
    for i in range(n):
        out_set = [j for j in range(n) if A[i][j] == 1]
        edges = sum(1 for j in out_set for k in out_set if A[j][k] == 1)
        profiles.append(edges)
    return tuple(sorted(profiles))

prof_T = out_edge_profile(A_nonnormal)
prof_Top = out_edge_profile(A_op)
print(f"  Out-edge profiles match: {prof_T == prof_Top}")
if prof_T != prof_Top:
    print(f"    T:    {prof_T}")
    print(f"    T^op: {prof_Top}")

# Try backtracking with smarter ordering
# Put most constrained vertex first: vertex with most edges to already-placed vertices
print("\n  Backtracking search with smart ordering...")

# BFS ordering from vertex 0
from collections import deque
visited = {0}
queue = deque([0])
bfs_order = [0]
while queue:
    v = queue.popleft()
    for u in range(n):
        if u not in visited and (A_nonnormal[v][u] == 1 or A_nonnormal[u][v] == 1):
            visited.add(u)
            queue.append(u)
            bfs_order.append(u)

result2 = backtrack_anti_auto(A_nonnormal, {}, bfs_order, n)
if result2 is not None:
    sigma = [result2[i] for i in range(n)]
    ok = all(A_nonnormal[sigma[i]][sigma[j]] == A_nonnormal[j][i]
             for i in range(n) for j in range(n) if i != j)
    print(f"  Anti-automorphism found! Verified: {ok}")
    print(f"  sigma = {sigma}")
else:
    print(f"  CONFIRMED: No anti-automorphism exists!")
    print(f"  This F_21 tournament is NOT self-converse!")

    # ADDITIONAL CHECK: Is T still vertex-transitive?
    # It should be, since it's a Cayley tournament.
    # Verify by checking that left multiplication by (1,0) is an automorphism.
    g0 = (1, 0)
    perm = [idx[mul(g0, elements[i])] for i in range(n)]
    aut_ok = all(A_nonnormal[perm[i]][perm[j]] == A_nonnormal[i][j]
                 for i in range(n) for j in range(n) if i != j)
    print(f"  Left mult by (1,0) is automorphism: {aut_ok}")

    g1 = (0, 1)
    perm1 = [idx[mul(g1, elements[i])] for i in range(n)]
    aut_ok1 = all(A_nonnormal[perm1[i]][perm1[j]] == A_nonnormal[i][j]
                  for i in range(n) for j in range(n) if i != j)
    print(f"  Left mult by (0,1) is automorphism: {aut_ok1}")

    print(f"\n  CONCLUSION: This is a vertex-transitive tournament that is NOT self-converse!")
    print(f"  THM-052's self-converse argument does NOT extend to all VT tournaments.")
    print(f"  We need a DIFFERENT proof strategy for non-abelian-group VT tournaments.")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
