#!/usr/bin/env python3
"""
CONNECTION: Irving-Omar det/per formula and our transfer matrix.

Irving-Omar (arXiv:2412.10572, Proposition 2):
  ham(D) = sum_{S⊆[n]} det(Ā[S]) · per(A[S^c])

where Ā = complement adjacency = A^T + I for tournaments.

Our transfer matrix:
  M[a,b] = sum_{S⊂U} (-1)^|S| · E_a(S∪{a}) · B_b(R∪{b})

QUESTION: Can we express M[a,b] using det/per of submatrices?

Observations:
1. E_a(W) ending at a = per of submatrix of A^T with endpoint constraint
2. B_b(W) starting at b = per of submatrix of A with endpoint constraint
3. The (-1)^|S| signs suggest det rather than per for one factor

HYPOTHESIS: M[a,b] is an endpoint-refinement of Irving-Omar's formula.

kind-pasteur-2026-03-06-S25 (continuation)
"""

from itertools import permutations
import numpy as np
from sympy import symbols, Matrix, Rational, expand, collect, Poly

def make_tournament(n, edges):
    """Create tournament from edge list. edges[(i,j)] = 1 means i->j."""
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = edges.get((i,j), 0)
    return T

def make_circulant(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def adjacency_matrix(T, n):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = T.get((i,j), 0)
    return A

def compute_M_entry(T, n, a, b):
    """Compute M[a,b] via definition."""
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


# ============================================================
# TEST: Irving-Omar det/per formula
# ============================================================
print("=" * 70)
print("Irving-Omar: ham(D) = sum_S det(Abar[S]) * per(A[S^c])")
print("=" * 70)

for n, name, S_set in [(3, "T_3", {1}), (5, "T_5 Paley", {1,2})]:
    T = make_circulant(n, S_set)
    A = adjacency_matrix(T, n)
    Abar = np.ones((n,n), dtype=int) - A  # Abar = J - A (includes self-loops!)

    # For tournament: Abar = A^T + I (since A[i,j]+A[j,i]=1 for i!=j, A[i,i]=0)
    assert np.array_equal(Abar, A.T + np.eye(n, dtype=int)), f"Abar != A^T + I at n={n}!"

    # Count ham paths directly
    H_direct = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(len(perm)-1):
            prod *= A[perm[k], perm[k+1]]
        H_direct += prod

    # Irving-Omar formula
    H_io = 0
    for mask in range(1 << n):
        S = [i for i in range(n) if mask & (1 << i)]
        Sc = [i for i in range(n) if not (mask & (1 << i))]
        if len(S) == 0:
            det_val = 1
        else:
            det_val = int(round(np.linalg.det(Abar[np.ix_(S, S)])))
        if len(Sc) == 0:
            per_val = 1
        else:
            # Compute permanent
            from itertools import permutations as perms
            per_val = 0
            for p in perms(range(len(Sc))):
                prod = 1
                for k in range(len(Sc)):
                    prod *= A[Sc[k], Sc[p[k]]]
                per_val += prod
        H_io += det_val * per_val

    print(f"\n  n={n} ({name}): H_direct = {H_direct}, H_IO = {H_io}, match: {H_direct == H_io}")
    print(f"    Ā = A^T: True")


# ============================================================
# KEY: Endpoint-conditioned det/per decomposition
# ============================================================
print("\n" + "=" * 70)
print("Endpoint-conditioned: H(a→b) via det/per?")
print("=" * 70)

# For a tournament, count paths from a to b
def count_paths_ab(A, n, a, b):
    """Count Hamiltonian paths from a to b."""
    count = 0
    for perm in permutations(range(n)):
        if perm[0] != a or perm[-1] != b:
            continue
        prod = 1
        for k in range(len(perm)-1):
            prod *= A[perm[k], perm[k+1]]
        count += prod
    return count

# Test: can we write H(a→b) as a det/per formula?
n = 5
T = make_circulant(n, {1,2})
A = adjacency_matrix(T, n)
Abar = A.T + np.eye(n, dtype=int)  # For tournament: Abar = A^T + I

print(f"\n  n={n}, Paley T_5:")
for a in range(min(n,3)):
    for b in range(min(n,3)):
        if a == b: continue
        hab = count_paths_ab(A, n, a, b)
        print(f"    H({a}->{b}) = {hab}", end="")

        # Try: H(a→b) = sum_S det(Ā[S∪{a}; rows/cols]) · per(A[Sc∪{b}; ...])
        # with endpoint constraints
        M_ab = compute_M_entry(T, n, a, b)
        print(f"  M[{a},{b}] = {M_ab}")

    print()


# ============================================================
# DISCOVERY: Relationship between M[a,b] and H(a→b)
# ============================================================
print("=" * 70)
print("Relationship between M[a,b] and H(a→b)")
print("=" * 70)

for n, name, S_set in [(3, "T_3", {1}), (5, "Paley_5", {1,2})]:
    T = make_circulant(n, S_set)
    A = adjacency_matrix(T, n)

    print(f"\n  n={n} ({name}):")

    M = np.zeros((n,n), dtype=int)
    Hab = np.zeros((n,n), dtype=int)
    for a in range(n):
        for b in range(n):
            M[a,b] = compute_M_entry(T, n, a, b)
            if a != b:
                Hab[a,b] = count_paths_ab(A, n, a, b)

    H = int(np.sum(Hab[0,:]))  # total from vertex 0

    print(f"    H(T) = {sum(Hab[i,j] for i in range(n) for j in range(n) if i != j)}")
    print(f"    M =")
    for row in M:
        print(f"      {list(row)}")
    print(f"    H(a→b) =")
    for row in Hab:
        print(f"      {list(row)}")
    print(f"    tr(M) = {int(np.trace(M))}")
    print(f"    sum(M) = {int(np.sum(M))}")

    # Check: is there a simple relation?
    # For VT tournaments: M = (H/n)*I, so M is NOT the H(a→b) matrix
    # But sum over endpoints: sum_a M[a,b] = column sum = sum_{a≠b} M[a,b] + M[b,b]
    # And sum_{a,b} M[a,b] = tr(M) + sum off-diag = H (for odd n)

    # Key question: what does M[a,b] for a≠b represent combinatorially?
    # Irving-Omar's det/per formula doesn't have endpoint conditioning
    # Our M[a,b] is the alternating sum of products of endpoint-conditioned counts

    # Check: M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
    # This is an INCLUSION-EXCLUSION over the internal vertices
    # It INCLUDES or EXCLUDES each internal vertex from the "left half" vs "right half"


# ============================================================
# CRITICAL TEST: Does sum_{a,b: a≠b} H(a→b) relate to Irving-Omar?
# ============================================================
print("\n" + "=" * 70)
print("Does M relate to the Irving-Omar det/per formula?")
print("=" * 70)

n = 5
T = make_circulant(n, {1,2})
A = adjacency_matrix(T, n)
Abar = A.T + np.eye(n, dtype=int)  # For tournament: Abar = A^T + I

print(f"\n  n={n} (Paley T_5):")
print(f"  A = ")
for row in A:
    print(f"    {list(row)}")

# Irving-Omar breaks ham(D) into S, S^c contributions
# Our M breaks it into (a, S, b) contributions
# Can we refine Irving-Omar to track endpoints?

# The permanent per(A[S^c]) counts CYCLE COVERS of S^c
# The determinant det(Ā[S]) counts SIGNED cycle covers of S using Ā = A^T

# A cycle cover is a collection of directed cycles covering all vertices
# A Hamiltonian PATH is NOT a cycle cover — it's a single path

# So Irving-Omar decomposes Ham path as: path = cycle cover decomposition via inclusion-exclusion
# The (-1)^|S| in our formula also comes from inclusion-exclusion

# INSIGHT: In the Irving-Omar formula, each term in the sum is:
#   det(A^T[S]) * per(A[S^c])
# = (sum_sigma sgn(sigma) prod A^T[i,sigma(i)] for i in S) * (sum_tau prod A[j,tau(j)] for j in S^c)
# = (sum_sigma sgn(sigma) prod A[sigma(i),i] for i in S) * (sum_tau prod A[j,tau(j)] for j in S^c)

# The first factor counts SIGNED cycle covers of S using REVERSE arcs
# The second factor counts cycle covers of S^c using forward arcs

# When S = {a} (single vertex), det(A^T[{a}]) = A^T[a,a] = 0 (no self-loops)
# So S with a single vertex contributes 0... unless we include the identity

# Wait, for tournaments A[i,i] = 0, so det(A^T[{i}]) = 0.
# But det(Ā[{i}]) = Ā[i,i] = 0 (since Ā = J-I-A, diagonal is 0)

# Hmm, the Irving-Omar formula includes S = empty set:
# det(Ā[∅]) = 1 (empty determinant), per(A[[n]]) = per(A) = #cycle covers
# So ham(D) = per(A) + sum_{S≠∅} det(Ā[S]) per(A[S^c])

# For tournaments: per(A) = # cycle covers of D
# And ham(D) = per(A) + correction terms from determinants

# Verify
from math import factorial

# Compute per(A) directly
def permanent(M):
    n = len(M)
    total = 0
    for p in permutations(range(n)):
        prod = 1
        for i in range(n):
            prod *= M[i][p[i]]
        total += prod
    return total

per_A = permanent(A)
print(f"  per(A) = {per_A}")

# H direct
H = 0
for perm in permutations(range(n)):
    prod = 1
    for k in range(len(perm)-1):
        prod *= A[perm[k], perm[k+1]]
    H += prod
print(f"  H(T) = {H}")

# Number of cycle covers
print(f"  Difference per(A) - H = {per_A - H}")
print(f"  (This is the contribution of non-Hamiltonian-path cycle covers)")


# ============================================================
# FORMULATE: Endpoint-conditioned Irving-Omar
# ============================================================
print("\n" + "=" * 70)
print("FORMULATION: M[a,b] as endpoint-conditioned Irving-Omar")
print("=" * 70)

print("""
The Irving-Omar formula decomposes ham(D) by splitting V into S and S^c:
  ham(D) = sum_S det(Ā[S]) · per(A[S^c])

Our transfer matrix decomposes by splitting V\\{a,b} into S and R:
  M[a,b] = sum_S (-1)^|S| · E_a(S+a) · B_b(R+b)

The ENDPOINT-CONDITIONED version of Irving-Omar would be:

  H(a→b) = sum_S det_*(Ā[S]; endpoint=a) · per_*(A[S^c]; endpoint=b)

where det_* and per_* are modified to track the endpoint structure.

KEY QUESTION: Can E_a(S+a) be expressed as a determinant or permanent
of a submatrix of Ā or A, with the endpoint constraint built in?

For E_a(W) = Ham paths through W ending at a:
  This is sum over permutations σ of W with σ(last) = a of product A[σ(i),σ(i+1)]
  = sum over orderings (w_1,...,w_m) with w_m = a of product_{k=1}^{m-1} A[w_k,w_{k+1}]

This is like a PERMANENT with row ordering constraint (last row = a).
In matrix terms: it's the (a, *) cofactor of a path-counting matrix.
""")


# ============================================================
# TEST: E_a(W) as modified permanent
# ============================================================
print("=" * 70)
print("TEST: E_a(W) via matrix permanent with endpoint constraint")
print("=" * 70)

n = 4
T = make_circulant(n, {1,2})
A = adjacency_matrix(T, n)

# For W = {0,1,2,3}, compute E_0(W) = paths ending at 0
W = list(range(n))
for a in range(n):
    ea = 0
    for perm in permutations(W):
        if perm[-1] != a: continue
        prod = 1
        for k in range(len(perm)-1):
            prod *= A[perm[k], perm[k+1]]
        ea += prod
    print(f"  E_{a}(V) = {ea}  (paths ending at {a})")

# Compare with column permanent or determinant structures
# The "last-vertex" constraint means we need specific matrix operations

print(f"\n  For comparison:")
print(f"    B_0(V) = {sum(1 for perm in permutations(W) if perm[0] == 0 and all(A[perm[k],perm[k+1]] for k in range(n-1)))}")

# At n=4 circulant {1,2}, the tournament is vertex-transitive
# So E_a = E_b for all a,b (and B_a = B_b) by VT symmetry


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
