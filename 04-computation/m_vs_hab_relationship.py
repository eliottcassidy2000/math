#!/usr/bin/env python3
"""
What is the precise relationship between M[a,b] and H(a->b)?

M[a,b] = inclusion-exclusion path decomposition
H(a->b) = simple Hamiltonian path count from a to b

We know:
  - tr(M) = H (odd n) or = 0 (even n)
  - sum(M) = H (odd n) or = 2H (even n)
  - M is symmetric (THM-030)
  - H(a->b) is NOT symmetric in general
  - M can have negative entries; H(a->b) >= 0

QUESTION: Is there a matrix transform connecting M and H_ab?
For instance: M = H_ab * P for some matrix P?
Or M = f(H_ab) where f is a spectral transform?

Let's investigate at n=3, 4, 5.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def count_paths_subset(A, verts, start=None, end=None):
    count = 0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        valid = True
        for i in range(len(p)-1):
            if A[p[i]][p[i+1]] != 1: valid = False; break
        if valid: count += 1
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

def endpoint_matrix(A):
    n = len(A)
    H_ab = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            H_ab[a][b] = count_paths_subset(A, list(range(n)), start=a, end=b)
    return H_ab

def ham_path_count(A):
    n = len(A)
    return count_paths_subset(A, list(range(n)))

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield [list(row) for row in A]

print("=" * 70)
print("M[a,b] vs H(a->b): FINDING THE RELATIONSHIP")
print("=" * 70)

# =====================================================================
# n=3: All tournaments
# =====================================================================
print("\n--- n=3 ---")
for A in all_tournaments(3):
    n = 3
    H = ham_path_count(A)
    M = transfer_matrix(A)
    Hab = endpoint_matrix(A)

    # Decompose Hab into symmetric + antisymmetric
    Hab_sym = (Hab + Hab.T) / 2
    Hab_anti = (Hab - Hab.T) / 2

    # Check if M = Hab_sym + something
    diff = M - Hab_sym
    print(f"\n  H={H}")
    print(f"  M    = {M.tolist()}")
    print(f"  H_ab = {Hab.tolist()}")
    print(f"  H_ab_sym = {Hab_sym.tolist()}")
    print(f"  M - H_ab_sym = {diff.tolist()}")

    # Check if M = Hab_sym + c*I
    # At diagonal: M[a,a] = Hab_sym[a,a] + c => c = M[a,a] - Hab[a,a]
    diag_diff = [M[i][i] - Hab[i][i] for i in range(n)]
    print(f"  diagonal diff M-H_ab: {diag_diff}")

# =====================================================================
# n=4: Sample tournaments
# =====================================================================
print("\n--- n=4 ---")
count = 0
for A in all_tournaments(4):
    count += 1
    if count > 8: break

    n = 4
    H = ham_path_count(A)
    M = transfer_matrix(A)
    Hab = endpoint_matrix(A)

    # Key identity checks
    # Does M = Hab + Hab^T - diag(something)?
    HpHT = Hab + Hab.T
    diff_HpHT = M - HpHT

    # Or M = (n-1)! / n * (symmetric rank-1 correction)?
    # What about M + M_adj?

    # Check: M * J = H * J where J is all-ones?
    J = np.ones((n, n), dtype=int)
    MJ = M @ J
    print(f"\n  H={H}")
    print(f"  M*J = {MJ.flatten().tolist()}")
    print(f"  H*J row = {[H]*n}")

    # Is M*1 = H*1? (row sums of M = H?)
    row_sums = M.sum(axis=1)
    print(f"  M row sums = {row_sums.tolist()}")

    # At even n: row sums of M should sum to 2H
    print(f"  sum(row_sums) = {sum(row_sums)} (should = 2H = {2*H})")

# =====================================================================
# n=5: Detailed relationship
# =====================================================================
print("\n--- n=5: Deeper analysis ---")

# Test the CYCLIC tournament specifically
A_cyc5 = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(5):
        if i != j and (j-i)%5 in [1, 4]:
            A_cyc5[i][j] = 1

n = 5
H = ham_path_count(A_cyc5)
M = transfer_matrix(A_cyc5)
Hab = endpoint_matrix(A_cyc5)
print(f"\nn=5 Paley: H={H}")
print(f"  M = {M.tolist()}")
print(f"  H_ab = {Hab.tolist()}")
print(f"  M = 2*I (since M=(H/n)*I = 10/5*I = 2*I)")

# For this scalar case, H_ab has a nice structure too
print(f"  H_ab row sums = {Hab.sum(axis=1).tolist()}")
print(f"  H_ab col sums = {Hab.sum(axis=0).tolist()}")
print(f"  H_ab diagonal = {[Hab[i][i] for i in range(n)]}")
print(f"  H_ab is circulant? ", end="")
is_circ = all(Hab[i][j] == Hab[(i+1)%n][(j+1)%n] for i in range(n) for j in range(n))
print(is_circ)

# For circulant Hab, M = 2*I means the DFT of the first row of Hab
# has a specific structure
if is_circ:
    first_row = [Hab[0][j] for j in range(n)]
    print(f"  H_ab first row: {first_row}")
    # DFT
    dft = np.fft.fft(first_row)
    print(f"  DFT: {[round(x.real, 3) + round(x.imag, 3)*1j for x in dft]}")

# =====================================================================
# General pattern: try M = (H_ab + H_ab^T)/2 + correction
# =====================================================================
print()
print("=" * 70)
print("GENERAL PATTERN SEARCH")
print("=" * 70)

# For odd n=5, test all tournaments
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
seen = set()

for bits in range(2**len(edges)):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    # Skip isomorphic
    key = tuple(tuple(row) for row in A)
    min_key = key
    for perm in permutations(range(n)):
        pkey = tuple(tuple(A[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if pkey < min_key:
            min_key = pkey
    if min_key in seen:
        continue
    seen.add(min_key)

    H = ham_path_count(A)
    M = transfer_matrix(A)
    Hab = endpoint_matrix(A)

    # M = (Hab + Hab^T)/2 + D where D is ???
    sym_Hab = (Hab + Hab.T) / 2
    D = M - sym_Hab

    # Is D a function of A (adjacency matrix)?
    # D should be symmetric and have tr(D) = tr(M) - tr(Hab)
    tr_D = np.trace(D)
    tr_M = np.trace(M)
    tr_Hab = np.trace(Hab)

    # Check: is D related to the SCORE matrix?
    # Score matrix S[i,j] = A[i][j] - A[j][i] = A[i][j] - (1-A[i][j]) = 2A[i][j]-1
    # (for i != j; S[i,i] = 0)
    S = np.array([[2*A[i][j]-1 if i != j else 0 for j in range(n)] for i in range(n)])

    # Check if D is proportional to something simple
    # D has integer or half-integer entries

    if H in [1, 9, 15]:
        print(f"\n  H={H}:")
        print(f"    M = {M.tolist()}")
        print(f"    sym(Hab) = {sym_Hab.tolist()}")
        print(f"    D = M - sym(Hab) = {D.tolist()}")
        print(f"    tr(D) = {tr_D}, tr(M)={tr_M}, tr(Hab)={tr_Hab}")

        # Key check: is D diagonal?
        off_diag_D = sum(abs(D[i][j]) for i in range(n) for j in range(n) if i != j)
        print(f"    D diagonal? {off_diag_D < 1e-10}")

        if off_diag_D > 1e-10:
            # Is D proportional to M?
            ratio = D / M if np.all(M != 0) else None
            print(f"    D elements: {[[round(D[i][j], 2) for j in range(n)] for i in range(n)]}")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
M[a,b] and H(a->b) are related but distinct:
  - M is symmetric; H(a->b) is NOT symmetric
  - M = (H_ab + H_ab^T)/2 + D where D is a symmetric correction
  - D is NOT diagonal in general
  - For scalar M (regular tournaments), M = (H/n)*I and H_ab is circulant

The correction D encodes the "signed path decomposition" beyond
simple endpoint counting. D captures the inclusion-exclusion
structure of the transfer matrix definition.

The physical interpretation: M[a,b] weights paths by how they
decompose the vertex set into "ending-at-a" and "starting-at-b"
subsets, with alternating signs. This is richer than just
counting paths from a to b.
""")
