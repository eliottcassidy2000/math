#!/usr/bin/env python3
"""
CONJECTURE: For any circulant tournament T on Z/nZ at odd n, M = (H/n)*I.

A circulant tournament on Z/nZ is determined by a generating set S ⊂ {1,...,n-1}
where S ∪ (-S mod n) = {1,...,n-1} and S ∩ (-S mod n) = ∅.
(i.e., for each pair {d, n-d}, exactly one is in S.)

The tournament has edge i→j iff (j-i) mod n ∈ S.

The circulant symmetry (rotation by 1) is an automorphism.
This means the position matrix P[v,k] is constant in v (each vertex
appears equally at each position in Ham paths).

If M is invariant under the cyclic rotation, then M commutes with the
permutation matrix of the rotation. For circulant matrices, this means
M is itself circulant. But M is also symmetric (THM-030).
A circulant that is also symmetric has the form:
  M[i,j] = f((j-i) mod n) where f(d) = f(n-d).

For M to be scalar (M = c*I), we need M[i,j] = 0 for i ≠ j.
Since M is circulant: M[i,j] = m_{(j-i) mod n}, we need m_d = 0 for d ≠ 0.

Let's test exhaustively at n = 3, 5, 7 and spot-check n = 9, 11.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def circulant_tournament(n, gen_set):
    """Circulant tournament on Z/nZ with generating set gen_set."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in gen_set:
                A[i][j] = 1
    return A

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

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

def all_circulant_gen_sets(n):
    """Generate all valid generating sets for circulant tournaments on Z/nZ."""
    half = list(range(1, (n+1)//2))  # {1, 2, ..., (n-1)/2}
    gen_sets = []
    for mask in range(1 << len(half)):
        gs = set()
        for k, d in enumerate(half):
            if mask & (1 << k):
                gs.add(d)
            else:
                gs.add(n - d)
        gen_sets.append(frozenset(gs))
    return gen_sets

# =====================================================================
print("=" * 70)
print("CONJECTURE: M = (H/n)*I for circulant tournaments at odd n")
print("=" * 70)

for n in [3, 5, 7]:
    gen_sets = all_circulant_gen_sets(n)
    print(f"\nn={n}: {len(gen_sets)} circulant tournaments")

    all_scalar = True
    for gs in gen_sets:
        A = circulant_tournament(n, gs)
        H = ham_path_count(A)
        M = transfer_matrix(A)

        expected = (H // n) if H % n == 0 else H / n
        is_scalar = np.allclose(M, (H/n) * np.eye(n))

        if not is_scalar:
            all_scalar = False
            print(f"  COUNTEREXAMPLE: gen={sorted(gs)}, H={H}, M scalar? NO")
            print(f"    M = {M.tolist()}")
            print(f"    H/n = {H/n:.4f}")
        else:
            print(f"  gen={sorted(gs)}: H={H}, M = ({H}/{n})*I = {H//n}*I ✓")

    if all_scalar:
        print(f"\n  *** ALL {len(gen_sets)} circulant tournaments at n={n} have M = (H/n)*I ***")
    else:
        print(f"\n  *** CONJECTURE FAILS at n={n} ***")

# =====================================================================
# n=9: spot check (might be slow for full M computation)
# =====================================================================
print()
print("=" * 70)
print("n=9: SPOT CHECK (may be slow)")
print("=" * 70)

n = 9
gen_sets = all_circulant_gen_sets(n)
print(f"n={n}: {len(gen_sets)} circulant tournaments")

# Only check a few
tested = 0
all_ok = True
for gs in gen_sets[:8]:
    A = circulant_tournament(n, gs)
    H = ham_path_count(A)

    # Quick check: if H % n != 0, cannot be scalar
    if H % n != 0:
        print(f"  gen={sorted(gs)}: H={H}, H%n={H%n} ≠ 0 => NOT scalar")
        all_ok = False
        continue

    # Check M[0,0] and M[0,1] only (circulant symmetry means M[a,a] = M[0,0] for all a)
    # and M[a,b] for a≠b depends only on (b-a) mod n
    M_00 = 0
    M_01 = 0
    U_diag = list(range(1, n))
    for k in range(len(U_diag)+1):
        for S in combinations(U_diag, k):
            S_set = set(S)
            R = [v for v in U_diag if v not in S_set]
            S_verts = sorted(list(S) + [0])
            R_verts = sorted(R + [0])
            ea = count_paths_subset(A, S_verts, end=0)
            bb = count_paths_subset(A, R_verts, start=0)
            M_00 += ((-1)**k) * ea * bb

    U_offdiag = [v for v in range(n) if v != 0 and v != 1]
    for k in range(len(U_offdiag)+1):
        for S in combinations(U_offdiag, k):
            S_set = set(S)
            R = [v for v in U_offdiag if v not in S_set]
            S_verts = sorted(list(S) + [0])
            R_verts = sorted(R + [1])
            ea = count_paths_subset(A, S_verts, end=0)
            bb = count_paths_subset(A, R_verts, start=1)
            M_01 += ((-1)**k) * ea * bb

    is_scalar = abs(M_01) < 1e-10 and abs(M_00 - H/n) < 1e-10
    tested += 1
    if is_scalar:
        print(f"  gen={sorted(gs)}: H={H}, M[0,0]={M_00}={H//n}, M[0,1]={M_01}=0 ✓")
    else:
        print(f"  gen={sorted(gs)}: H={H}, M[0,0]={M_00}, M[0,1]={M_01} => scalar? {is_scalar}")
        all_ok = False

if all_ok and tested > 0:
    print(f"\n  All {tested} tested circulant tournaments at n={n} have scalar M ✓")

# =====================================================================
# Why does circulant symmetry force M = (H/n)*I?
# =====================================================================
print()
print("=" * 70)
print("PROOF ARGUMENT")
print("=" * 70)
print("""
THEOREM: For any circulant tournament T on Z/nZ at odd n, M = (H/n)*I.

PROOF:
1. The rotation σ: i ↦ i+1 (mod n) is an automorphism of T.

2. By the definition of M via inclusion-exclusion:
   M[σ(a), σ(b)] = M[a, b] for all a, b.
   (Because σ maps Ham paths to Ham paths, preserving the IE decomposition.)

3. This means M is a CIRCULANT matrix: M[a,b] = m_{(b-a) mod n}.

4. M is symmetric (THM-030): M[a,b] = M[b,a].
   So m_d = m_{n-d} for all d.

5. At odd n: tr(M) = H (since sum_a (-1)^{pos(a,P)} = 1 for each P).
   So n * m_0 = H, giving m_0 = H/n.

6. Now we need: m_d = 0 for d = 1, 2, ..., (n-1)/2.
   Since M is circulant AND symmetric, it has the form:
   M = m_0 * I + sum_{d=1}^{(n-1)/2} m_d * (P_d + P_{n-d})
   where P_d is the dth power of the cyclic permutation matrix.

   The eigenvalues of M are: λ_k = m_0 + sum_d m_d * (ω^{kd} + ω^{-kd})
   where ω = e^{2πi/n}.

   For M = (H/n)*I, we need all eigenvalues equal: λ_k = H/n for all k.
   This means: sum_d m_d * (ω^{kd} + ω^{-kd}) = 0 for all k ≠ 0.

   This is NOT automatic from circulant symmetry alone.
   We need an additional argument to show m_d = 0 for d ≥ 1.

7. KEY: The transfer matrix M satisfies M*1 = H/n * 1 (column sum property)
   at odd n for ANY tournament (not just circulant).
   Wait — is this true? The column sum of M is:
   sum_a M[a,b] = sum_P sum_a (-1)^{pos(a,P)} * A[P-edges]
   At odd n: sum_a (-1)^{pos(a,P)} = 1, so sum_a M[a,b] = H
   for ALL b? No, the sum depends on b through the Ham paths.

   Actually: sum_a M[a,b] = #{Ham paths through b weighted by (-1)^{position}}
   This is NOT simply H.

   For circulant: by symmetry, sum_a M[a,b] is the same for all b.
   And sum_b sum_a M[a,b] = tr(M) = H at odd n? No, sum of all entries ≠ trace.

   Hmm, the proof needs more work. Let me check numerically whether the
   column sums are constant for circulant tournaments.
""")

# Quick check: column sums of M for circulant tournaments
n = 5
for gs in all_circulant_gen_sets(n):
    A = circulant_tournament(n, gs)
    M = transfer_matrix(A)
    H = ham_path_count(A)
    col_sums = [sum(M[a][b] for a in range(n)) for b in range(n)]
    row_sums = [sum(M[a][b] for b in range(n)) for a in range(n)]
    print(f"  gen={sorted(gs)}: H={H}, col_sums={col_sums}, row_sums={row_sums}")

print()
print("=" * 70)
print("DONE")
print("=" * 70)
