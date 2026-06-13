#!/usr/bin/env python3
"""
When is M = (H/n)*I? (Scalar transfer matrix)

At n=3: the cyclic tournament (regular) has M = I = (3/3)*I.
At n=5: BOTH class 9 (score [3,2,2,2,1]) and class 11 (score [2,2,2,2,2])
have M = 3*I = (15/5)*I.

Class 11 is the Paley/regular tournament — expected by vertex-transitivity.
Class 9 is NOT regular — what property makes it scalar?

HYPOTHESIS: M = (H/n)*I iff T is "path-transitive" in some sense:
every vertex plays the same role in Hamiltonian path counting.

This means M[a,a] = H/n for all a, i.e., the "diagonal path count"
is uniform. And M[a,b] = 0 for all a != b.

M[a,a] = sum_S (-1)^|S| E_a(S+a) B_a(R+a) counts something like
"how special is vertex a in the path structure."
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start:
            continue
        if end is not None and p[-1] != end:
            continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

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

def ham_path_count(A):
    n = len(A)
    return count_paths_subset(A, list(range(n)))

def tournament_canonical(A):
    n = len(A)
    min_adj = None
    for perm in permutations(range(n)):
        adj = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if min_adj is None or adj < min_adj:
            min_adj = adj
    return min_adj

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield bits, [list(row) for row in A]

def automorphism_group_size(A):
    """Count automorphisms of tournament A."""
    n = len(A)
    count = 0
    for perm in permutations(range(n)):
        is_auto = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != A[perm[i]][perm[j]]:
                    is_auto = False
                    break
            if not is_auto:
                break
        if is_auto:
            count += 1
    return count

# =====================================================================
# n=5: Deep dive into the scalar classes
# =====================================================================
print("=" * 70)
print("SCALAR M AT n=5: Which tournaments have M = (H/n)*I?")
print("=" * 70)

n = 5
seen_canonical = {}
scalar_classes = []
nonscalar_classes = []

for bits, A in all_tournaments(n):
    canon = tournament_canonical(A)
    if canon in seen_canonical:
        continue

    H = ham_path_count(A)
    M = transfer_matrix(A)
    is_scalar = np.allclose(M, (H/n)*np.eye(n))

    scores = sorted([sum(row) for row in A], reverse=True)
    auto_size = automorphism_group_size(A)

    seen_canonical[canon] = {
        'A': A, 'H': H, 'M': M, 'scores': scores,
        'is_scalar': is_scalar, 'auto_size': auto_size
    }

    if is_scalar:
        scalar_classes.append(canon)
    else:
        nonscalar_classes.append(canon)

print(f"\nTotal isomorphism classes at n={n}: {len(seen_canonical)}")
print(f"Scalar (M = (H/n)*I): {len(scalar_classes)}")
print(f"Non-scalar: {len(nonscalar_classes)}")

print(f"\nScalar classes:")
for canon in scalar_classes:
    d = seen_canonical[canon]
    print(f"  scores={d['scores']}, H={d['H']}, |Aut|={d['auto_size']}, H/n={d['H']/n}")

print(f"\nNon-scalar classes (showing diagonal of M):")
for canon in sorted(nonscalar_classes, key=lambda c: seen_canonical[c]['H']):
    d = seen_canonical[canon]
    diag = [d['M'][i][i] for i in range(n)]
    print(f"  scores={d['scores']}, H={d['H']}, |Aut|={d['auto_size']}, diag={diag}")

# =====================================================================
# Property analysis: what do scalar classes share?
# =====================================================================
print()
print("=" * 70)
print("PROPERTY ANALYSIS OF SCALAR CLASSES")
print("=" * 70)

for canon in scalar_classes:
    d = seen_canonical[canon]
    A = d['A']
    n_loc = len(A)

    # Source/sink count
    outdeg = [sum(row) for row in A]
    indeg = [n_loc - 1 - od for od in outdeg]

    # Endpoint path counts
    start_counts = [count_paths_subset(A, list(range(n_loc)), start=v) for v in range(n_loc)]
    end_counts = [count_paths_subset(A, list(range(n_loc)), end=v) for v in range(n_loc)]

    print(f"\n  scores={d['scores']}, H={d['H']}, |Aut|={d['auto_size']}")
    print(f"  outdegrees = {outdeg}")
    print(f"  #paths starting at v = {start_counts}")
    print(f"  #paths ending at v = {end_counts}")
    print(f"  start counts uniform? {len(set(start_counts)) == 1}")
    print(f"  end counts uniform? {len(set(end_counts)) == 1}")

    # Check: is every vertex in the same number of 3-cycles?
    tri_count = [0] * n_loc
    for i in range(n_loc):
        for j in range(n_loc):
            for k in range(n_loc):
                if i != j and j != k and i != k:
                    if A[i][j] == 1 and A[j][k] == 1 and A[k][i] == 1:
                        tri_count[i] += 1
    tri_per_vertex = [t // 2 for t in tri_count]  # each 3-cycle counted twice per vertex
    print(f"  3-cycles per vertex = {tri_per_vertex}")

# =====================================================================
# n=3: Verify (should be only cyclic = regular)
# =====================================================================
print()
print("=" * 70)
print("SCALAR M AT n=3")
print("=" * 70)

n = 3
seen3 = {}
for bits, A in all_tournaments(n):
    canon = tournament_canonical(A)
    if canon in seen3:
        continue
    H = ham_path_count(A)
    M = transfer_matrix(A)
    is_scalar = np.allclose(M, (H/n)*np.eye(n))
    scores = sorted([sum(row) for row in A], reverse=True)
    auto_size = automorphism_group_size(A)
    seen3[canon] = True
    label = "SCALAR" if is_scalar else ""
    m_str = f"{int(H/n)}*I" if is_scalar else str(M.tolist())
    print(f"  scores={scores}, H={H}, |Aut|={auto_size}, M={m_str} {label}")

# =====================================================================
# n=7: Check which classes are scalar (expensive!)
# =====================================================================
print()
print("=" * 70)
print("SCALAR M AT n=7 (sampling)")
print("=" * 70)

# n=7 has C(7,2)=21 edges, so 2^21 = 2M tournaments.
# Too many to enumerate, but we can sample and check specific ones.

n = 7

# Paley tournament: QR mod 7 = {1, 2, 4}
A_paley = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if i != j:
            if (j - i) % 7 in [1, 2, 4]:
                A_paley[i][j] = 1

H_paley = ham_path_count(A_paley)
M_paley = transfer_matrix(A_paley)
is_scalar_paley = np.allclose(M_paley, (H_paley/7)*np.eye(7))
scores_paley = sorted([sum(row) for row in A_paley], reverse=True)
print(f"\n  Paley n=7: scores={scores_paley}, H={H_paley}, scalar={is_scalar_paley}")
if is_scalar_paley:
    print(f"  M = {H_paley/7:.0f} * I")

# Circulant tournament: i->j if (j-i) mod 7 in {1, 2, 3}
A_circ = [[0]*7 for _ in range(7)]
for i in range(7):
    for j in range(7):
        if i != j:
            if (j - i) % 7 in [1, 2, 3]:
                A_circ[i][j] = 1

H_circ = ham_path_count(A_circ)
M_circ = transfer_matrix(A_circ)
is_scalar_circ = np.allclose(M_circ, (H_circ/7)*np.eye(7))
scores_circ = sorted([sum(row) for row in A_circ], reverse=True)
print(f"\n  Circulant {1,2,3} n=7: scores={scores_circ}, H={H_circ}, scalar={is_scalar_circ}")
if is_scalar_circ:
    print(f"  M = {H_circ/7:.0f} * I")
else:
    diag = [M_circ[i][i] for i in range(7)]
    print(f"  diagonal = {diag}")
    print(f"  M[0,1] = {M_circ[0][1]}")

# Near-regular: one vertex with different score
# Take Paley and flip one arc
A_near = [row[:] for row in A_paley]
A_near[0][1] = 0  # was 1
A_near[1][0] = 1  # was 0
H_near = ham_path_count(A_near)
M_near = transfer_matrix(A_near)
is_scalar_near = np.allclose(M_near, (H_near/7)*np.eye(7))
scores_near = sorted([sum(row) for row in A_near], reverse=True)
print(f"\n  Near-Paley n=7: scores={scores_near}, H={H_near}, scalar={is_scalar_near}")
if not is_scalar_near:
    diag = [M_near[i][i] for i in range(7)]
    print(f"  diagonal = {diag}")

print()
print("=" * 70)
print("CONJECTURE")
print("=" * 70)
print("""
M = (H/n)*I seems to hold when EVERY vertex plays the same role
in Hamiltonian path counting:
  - same number of paths starting/ending at each vertex
  - same "path centrality" in the inclusion-exclusion

At n=5, this happens for:
  1. The regular tournament (score [2,2,2,2,2]) — vertex-transitive
  2. A non-regular tournament (score [3,2,2,2,1]) — NOT vertex-transitive
     but still "path-symmetric"

This second case is surprising! It means vertex-transitivity is
SUFFICIENT but NOT NECESSARY for M to be scalar.

QUESTION: Is there a weaker property than vertex-transitivity that
characterizes when M is scalar? Perhaps "path-transitivity" where
the automorphism group acts transitively on Hamiltonian paths?
""")
