"""
vitali_dc5_algebraic_proof.py -- kind-pasteur-2026-03-13-S61

Attempt to prove dc5 = 0 algebraically.

Key observation: dc5 = 0 is about the TOTAL directed 5-cycle count,
summed over ALL vertex sets. This is a global invariant.

Approach: Express dc5 as a LINEAR FUNCTION of the arc changes within S.

For a tournament T, define:
  f_5(T) = total number of directed 5-cycles in T
         = sum over all 5-vertex V of hc(T[V])

The reversal on S changes arcs a_{ij} -> 1-a_{ij} for i,j in S (i != j).
These are 6 arcs (12 ordered pairs, but 6 independent bits).

CLAIM: f_5(T') - f_5(T) is a LINEAR function of the arc changes.

This is because each directed 5-cycle is a PRODUCT of 5 arc indicators.
The change in a product when one factor changes is:
  prod_after - prod_before = (factor change) * (product of other factors)

But if MULTIPLE factors change simultaneously (for k=3,4 arcs in S),
the difference is NOT linear — it's polynomial of degree k.

So the linearity claim is FALSE in general. But for |V cap S| = 2,
only 1 arc within S is in V, so the change IS linear.

For |V cap S| = 3,4, the change is polynomial of degree 3 or 6.

But dc5 = 0 regardless of the external structure. This suggests
a CANCELLATION at the level of the (1,1,2,2) complement.

Let me try a different approach: COUNT the total directed cycles
as a polynomial in the 6 S-arcs, keeping all other arcs fixed.

f_5(T) as a function of a_{ij} for (i,j) in S x S is a multilinear
polynomial (since each product of arc indicators is multilinear).
The complement maps a_{ij} -> 1 - a_{ij}. Under this transformation,
a multilinear function f transforms as:
  f(1-x) = ... (expansion)

For a linear function f(x) = c_0 + c_1 * x:
  f(1-x) = c_0 + c_1 - c_1 * x

So f(1-x) - f(x) = c_1 - 2*c_1*x = c_1*(1-2x).
If the (1,1,2,2) score constraint means that the ARC VALUES at the
specific (1,1,2,2) tournament satisfy certain identities, this
could force the sum to vanish.

Actually, let's be concrete. The 4-vertex tournament T[S] on S has
3 distinct (1,1,2,2) tournaments up to labeling. But there's exactly
ONE tournament on 4 vertices with score (1,1,2,2) up to isomorphism.

Wait, actually on LABELED vertices {s1,s2,s3,s4}, there are multiple
labelings that give score (1,1,2,2). Let me enumerate:

For 4 vertices, out-degrees (d1,d2,d3,d4) is a permutation of (1,1,2,2).
The number of tournaments with a GIVEN degree sequence is:
C(4,2)!/prod... actually this is complicated.

Let me instead verify the dc5=0 identity by expressing dc5 as a
POLYNOMIAL in the S-arcs and checking that it vanishes at all
(1,1,2,2) points and their complements.

We have 6 S-arc variables: a_{01}, a_{02}, a_{03}, a_{12}, a_{13}, a_{23}
(where S = {s0, s1, s2, s3} and a_{ij} = 1 iff s_i -> s_j).

The constraint is: the TOTAL out-degree of s_i within S is d_i,
where (d0,d1,d2,d3) is a permutation of (1,1,2,2).

The complement sends a_{ij} -> 1-a_{ij}, changing degrees to (3-d_i).

dc5 = f_5(complement) - f_5(original) = sum over all 5-vertex V of
  [hc(T'[V]) - hc(T[V])]

For a FIXED external configuration (arcs outside S), dc5 is a function
of the 6 S-arcs only. We need to show it's 0 at all (1,1,2,2) points.

Let me compute this for a specific small external configuration.
Take n=7 (smallest n where 5-cycles exist and |V cap S| can be 2).
Wait, at n=7 we have |S|=4 and 3 external vertices.
5-vertex sets:
  |V cap S| = 2: C(4,2)*C(3,3) = 6 sets (each has 2 S-verts + all 3 ext)
  |V cap S| = 3: C(4,3)*C(3,2) = 12 sets
  |V cap S| = 4: C(4,4)*C(3,1) = 3 sets

For each set V, hc(T[V]) is a polynomial in the S-arcs.

This is a finite computation. Let me do it symbolically.
"""

import numpy as np
from itertools import combinations, permutations, product
from collections import Counter

def count_directed_5cycles_on_set(A, vset):
    """Count directed 5-cycles on vertex set (fixing first vertex)."""
    combo = tuple(sorted(vset))
    count = 0
    for perm in permutations(combo[1:]):
        path = (combo[0],) + perm
        valid = True
        for i in range(5):
            if A[path[i]][path[(i+1) % 5]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

# Enumerate all (1,1,2,2) tournaments on 4 vertices
# S = {0,1,2,3}, 6 arcs: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
print("All (1,1,2,2) labeled tournaments on {0,1,2,3}:")
print("=" * 60)

labeled_1122 = []
for bits in range(64):
    A = np.zeros((4, 4), dtype=int)
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    scores = tuple(sorted([sum(A[i]) for i in range(4)]))
    if scores == (1, 1, 2, 2):
        labeled_1122.append(A.copy())
        score_actual = tuple(sum(A[i]) for i in range(4))
        arcs = []
        for i in range(4):
            for j in range(i+1, 4):
                if A[i][j]:
                    arcs.append(f"{i}->{j}")
                else:
                    arcs.append(f"{j}->{i}")
        print(f"  scores={score_actual}, arcs={arcs}")

print(f"\n  Total: {len(labeled_1122)} labeled (1,1,2,2) tournaments")

# Verify complement pairing
print(f"\n  Complement pairing:")
for idx, A in enumerate(labeled_1122):
    # Complement: flip all arcs
    B = 1 - A
    for i in range(4):
        B[i][i] = 0
    scores_B = tuple(sorted([sum(B[i]) for i in range(4)]))
    # Find B in the list
    found = -1
    for jdx, C in enumerate(labeled_1122):
        if np.array_equal(B, C):
            found = jdx
            break
    print(f"  T[{idx}] complement = T[{found}]")

# Now: for n=7 with S={0,1,2,3} and ext={4,5,6},
# fix the external arcs and compute dc5 as function of S-arcs.
n = 7
S = [0, 1, 2, 3]
ext = [4, 5, 6]

print(f"\n{'='*60}")
print(f"dc5 AS FUNCTION OF S-ARCS (n={n})")
print(f"{'='*60}")

# Fix a specific external configuration
# ext arcs: 4->5, 4->6, 5->6 (transitive on ext) and
# cross arcs: fix some configuration
np.random.seed(42)

for config in range(5):
    # Random cross arcs (S <-> ext)
    cross = np.zeros((n, n), dtype=int)
    for s in S:
        for e in ext:
            if np.random.random() < 0.5:
                cross[s][e] = 1
            else:
                cross[e][s] = 1

    # Random ext-ext arcs
    for i in range(len(ext)):
        for j in range(i+1, len(ext)):
            if np.random.random() < 0.5:
                cross[ext[i]][ext[j]] = 1
            else:
                cross[ext[j]][ext[i]] = 1

    # Now compute dc5 for each (1,1,2,2) S-tournament and its complement
    results = []
    for T_S in labeled_1122:
        A = cross.copy()
        for i in range(4):
            for j in range(4):
                if i != j:
                    A[S[i]][S[j]] = T_S[i][j]

        # Count total directed 5-cycles
        total_5 = 0
        for combo in combinations(range(n), 5):
            total_5 += count_directed_5cycles_on_set(A, combo)
        results.append(total_5)

    # Check: for each (T, complement(T)) pair, is dc5 = 0?
    print(f"\n  Config {config}:")
    for idx in range(len(labeled_1122)):
        # Find complement
        B = 1 - labeled_1122[idx]
        for i in range(4):
            B[i][i] = 0
        comp_idx = -1
        for jdx, C in enumerate(labeled_1122):
            if np.array_equal(B, C):
                comp_idx = jdx
                break
        dc5 = results[comp_idx] - results[idx]

        if dc5 != 0:
            print(f"    T[{idx}] -> T[{comp_idx}]: f5={results[idx]} -> {results[comp_idx]}, dc5={dc5}")
        else:
            pass  # Don't print zeros to keep output clean

    # But we also need lambda-preservation. Check which pairs have it.
    print(f"  Lambda-preserving pairs:")
    for idx in range(len(labeled_1122)):
        B_S = 1 - labeled_1122[idx]
        for i in range(4):
            B_S[i][i] = 0
        comp_idx = -1
        for jdx, C in enumerate(labeled_1122):
            if np.array_equal(B_S, C):
                comp_idx = jdx
                break

        # Build full A and B
        A = cross.copy()
        B = cross.copy()
        for i in range(4):
            for j in range(4):
                if i != j:
                    A[S[i]][S[j]] = labeled_1122[idx][i][j]
                    B[S[i]][S[j]] = labeled_1122[comp_idx][i][j]

        # Check lambda
        lam_A = np.zeros((n, n), dtype=int)
        lam_B = np.zeros((n, n), dtype=int)
        for u in range(n):
            for v in range(u+1, n):
                for w in range(n):
                    if w == u or w == v:
                        continue
                    if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                        lam_A[u][v] += 1
                    if (B[u][v] and B[v][w] and B[w][u]) or (B[v][u] and B[u][w] and B[w][v]):
                        lam_B[u][v] += 1

        lam_preserved = np.array_equal(lam_A, lam_B)

        f5_A = results[idx]
        f5_B = results[comp_idx]
        dc5 = f5_B - f5_A

        if lam_preserved:
            print(f"    T[{idx}] -> T[{comp_idx}]: dc5={dc5}, lambda PRESERVED")

# Key check: does dc5=0 ONLY when lambda is preserved, or always?
print(f"\n{'='*60}")
print(f"dc5 WITHOUT lambda preservation?")
print(f"{'='*60}")

np.random.seed(999)
dc5_with_lambda = []
dc5_without_lambda = []

for config in range(100):
    cross = np.zeros((n, n), dtype=int)
    for s in S:
        for e in ext:
            cross[s][e] = int(np.random.random() < 0.5)
            cross[e][s] = 1 - cross[s][e]
    for i in range(len(ext)):
        for j in range(i+1, len(ext)):
            cross[ext[i]][ext[j]] = int(np.random.random() < 0.5)
            cross[ext[j]][ext[i]] = 1 - cross[ext[i]][ext[j]]

    results = []
    for T_S in labeled_1122:
        A = cross.copy()
        for i in range(4):
            for j in range(4):
                if i != j:
                    A[S[i]][S[j]] = T_S[i][j]
        total_5 = sum(count_directed_5cycles_on_set(A, combo) for combo in combinations(range(n), 5))
        results.append(total_5)

    for idx in range(len(labeled_1122)):
        B_S = 1 - labeled_1122[idx]
        for i in range(4):
            B_S[i][i] = 0
        comp_idx = next(j for j, C in enumerate(labeled_1122) if np.array_equal(B_S, C))

        A_full = cross.copy()
        B_full = cross.copy()
        for i in range(4):
            for j in range(4):
                if i != j:
                    A_full[S[i]][S[j]] = labeled_1122[idx][i][j]
                    B_full[S[i]][S[j]] = labeled_1122[comp_idx][i][j]

        lam_A = np.zeros((n, n), dtype=int)
        lam_B = np.zeros((n, n), dtype=int)
        for u in range(n):
            for v in range(u+1, n):
                for w in range(n):
                    if w == u or w == v:
                        continue
                    if (A_full[u][v] and A_full[v][w] and A_full[w][u]) or (A_full[v][u] and A_full[u][w] and A_full[w][v]):
                        lam_A[u][v] += 1
                    if (B_full[u][v] and B_full[v][w] and B_full[w][u]) or (B_full[v][u] and B_full[u][w] and B_full[w][v]):
                        lam_B[u][v] += 1

        lam_preserved = np.array_equal(lam_A, lam_B)
        dc5 = results[comp_idx] - results[idx]

        if lam_preserved:
            dc5_with_lambda.append(dc5)
        else:
            dc5_without_lambda.append(dc5)

print(f"  Lambda-preserved: {len(dc5_with_lambda)} reversals")
print(f"    dc5 distribution: {dict(sorted(Counter(dc5_with_lambda).items()))}")
print(f"    dc5 = 0 always? {all(x == 0 for x in dc5_with_lambda)}")

print(f"\n  Lambda NOT preserved: {len(dc5_without_lambda)} reversals")
print(f"    dc5 distribution: {dict(sorted(Counter(dc5_without_lambda).items()))}")
print(f"    dc5 = 0 always? {all(x == 0 for x in dc5_without_lambda)}")

print("\nDone.")
