#!/usr/bin/env python3
"""
Spectral analysis of the symmetric transfer matrix M.

Since M[a,b] = M[b,a] (THM-030), M is a real symmetric matrix,
so all eigenvalues are real. What is the spectral structure?

Key questions:
1. What are the eigenvalues for specific tournaments?
2. Is there a relationship between eigenvalues and H(T)?
3. Does the all-ones vector have a special role (row sums)?
4. What is the characteristic polynomial of M?
"""

import numpy as np
from itertools import permutations, combinations
from fractions import Fraction

def tournament_from_matrix(A):
    """A is the adjacency matrix (A[i][j] = 1 if i->j)."""
    n = len(A)
    return A

def ham_paths_count(A, start=None, end=None):
    """Count Hamiltonian paths in tournament A."""
    n = len(A)
    count = 0
    for p in permutations(range(n)):
        if start is not None and p[0] != start:
            continue
        if end is not None and p[-1] != end:
            continue
        valid = True
        for i in range(n-1):
            if A[p[i]][p[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def transfer_matrix(A):
    """Compute transfer matrix M[a,b] for tournament A."""
    n = len(A)
    M = np.zeros((n, n))
    vertices = list(range(n))

    for a in range(n):
        for b in range(n):
            U = [v for v in vertices if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]

                    # E_a(S + {a}): paths ending at a through S + {a}
                    S_verts = sorted(list(S) + [a])
                    ea = ham_paths_count(A, end=a) if len(S_verts) == n else count_paths_through(A, S_verts, end=a)

                    # B_b(R + {b}): paths starting at b through R + {b}
                    R_verts = sorted(R + [b])
                    bb = count_paths_through(A, R_verts, start=b)

                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def count_paths_through(A, verts, start=None, end=None):
    """Count Hamiltonian paths through exactly the vertices in verts."""
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

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(edges)):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield [list(row) for row in A]

print("=" * 70)
print("SPECTRAL ANALYSIS OF TRANSFER MATRIX M")
print("=" * 70)

# n=3: two non-isomorphic tournaments
print("\n--- n=3 ---")
# Transitive: 0->1->2, 0->2
A_trans3 = [[0,1,1],[0,0,1],[0,0,0]]
# Cyclic: 0->1->2->0
A_cyc3 = [[0,1,0],[0,0,1],[1,0,0]]

for name, A in [("Transitive", A_trans3), ("3-cycle", A_cyc3)]:
    M = transfer_matrix(A)
    H = ham_paths_count(A)
    evals = np.linalg.eigvalsh(M)
    print(f"  {name}: H={H}, M = {M.astype(int).tolist()}")
    print(f"    eigenvalues: {sorted(evals)}")
    print(f"    trace = {np.trace(M):.0f}, sum = {M.sum():.0f}")
    print(f"    row sums = {M.sum(axis=1).astype(int).tolist()}")

# n=4
print("\n--- n=4 ---")
# Collect unique (H, spectrum) pairs
from collections import defaultdict
h_spectra = defaultdict(list)

for A in all_tournaments(4):
    M = transfer_matrix(A)
    H = ham_paths_count(A)
    evals = sorted(np.linalg.eigvalsh(M))
    evals_rounded = tuple(round(e, 6) for e in evals)
    h_spectra[(H, evals_rounded)] = M

print(f"  Distinct (H, spectrum) pairs: {len(h_spectra)}")
for (H, evals), M in sorted(h_spectra.items()):
    print(f"    H={H}: eigenvalues = {list(evals)}")
    print(f"      trace={np.trace(M):.0f}, off-diag sum={M.sum()-np.trace(M):.0f}")
    print(f"      row sums = {M.sum(axis=1).astype(int).tolist()}")

# n=5
print("\n--- n=5 (sample) ---")
count_by_H = defaultdict(int)
spectra_by_H = defaultdict(set)

for A in all_tournaments(5):
    M = transfer_matrix(A)
    H = ham_paths_count(A)
    evals = sorted(np.linalg.eigvalsh(M))
    evals_rounded = tuple(round(e, 4) for e in evals)
    count_by_H[H] += 1
    spectra_by_H[H].add(evals_rounded)

print(f"  H values: {sorted(count_by_H.keys())}")
for H in sorted(spectra_by_H.keys()):
    specs = spectra_by_H[H]
    print(f"  H={H} ({count_by_H[H]} tournaments, {len(specs)} distinct spectra):")
    for s in sorted(specs)[:3]:  # show up to 3
        print(f"    {list(s)}")

# Key insight check: is the largest eigenvalue related to H?
print("\n--- Eigenvalue-H relationship ---")
print("  For n=4:")
for (H, evals), M in sorted(h_spectra.items()):
    max_eval = max(evals)
    trace = sum(evals)
    total = M.sum()
    print(f"    H={H}: max_eval={max_eval}, trace={trace}, total_sum={total}")

# Check: for odd n, T = tr(M) = H and Sigma = 0
# For odd n, row sums should all be equal if Sigma = 0?
# No — Sigma = 0 means SUM of off-diag = 0, not that each row sum is the same.
print("\n--- Row sum structure ---")
print("  For odd n=3:")
for name, A in [("Transitive", A_trans3), ("3-cycle", A_cyc3)]:
    M = transfer_matrix(A)
    rs = M.sum(axis=1)
    print(f"  {name}: row sums = {rs}")
    print(f"    M[a,a] = {[M[a,a] for a in range(3)]}")

print("\n  For odd n=5 (first few):")
count = 0
for A in all_tournaments(5):
    if count >= 5:
        break
    M = transfer_matrix(A)
    H = ham_paths_count(A)
    rs = M.sum(axis=1)
    diag = [M[a][a] for a in range(5)]
    print(f"    H={H}: row_sums={[round(x,1) for x in rs]}, diag={[round(x,1) for x in diag]}")
    count += 1

# For even n=4:
print("\n  For even n=4 (first few):")
count = 0
for A in all_tournaments(4):
    if count >= 5:
        break
    M = transfer_matrix(A)
    H = ham_paths_count(A)
    rs = M.sum(axis=1)
    diag = [M[a][a] for a in range(4)]
    off_diag_sum = M.sum() - np.trace(M)
    print(f"    H={H}: row_sums={[round(x,1) for x in rs]}, diag={diag}, Sigma={off_diag_sum:.0f}")
    count += 1

# The KEY relationship for even n: H = (1/2)*Sigma (at c=1, r=1/2)
# And tr(M) = 0
print("\n--- Even n: H = Sigma/2, tr(M) = 0 ---")
for A in all_tournaments(4):
    M = transfer_matrix(A)
    H = ham_paths_count(A)
    Sigma = M.sum() - np.trace(M)
    tr = np.trace(M)
    assert abs(tr) < 1e-10, f"tr(M) = {tr} ≠ 0 for even n"
    assert abs(H - Sigma/2) < 1e-10, f"H = {H} ≠ Sigma/2 = {Sigma/2}"
print("  ALL n=4 tournaments: tr(M)=0 and H=Sigma/2. VERIFIED.")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
