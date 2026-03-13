#!/usr/bin/env python3
"""
paley_poincare_duality.py — opus-2026-03-13-S71

Investigate Poincaré duality for Paley tournament path complexes.

Key observation: For Paley P_p (p ≡ 3 mod 4):
  Ω_m/p is palindromic around m = (p-2)/2

This suggests the GLMY path complex of P_p has a duality structure.

For a complete symmetric graph K_n (not tournament), the clique complex
has Poincaré duality. For a tournament, we need a replacement concept.

Key idea: The Paley tournament on Z_p has the property that
  QR ∪ NQR = Z_p^* (complete, no overlap for p ≡ 3 mod 4)
  The QR graph is self-complementary: P_p ≅ P_p^op (up to a permutation!)
  Specifically, multiplication by a non-residue is an anti-automorphism.

Could the self-complementarity of P_p give rise to Poincaré duality?

Also: does the PATH REVERSAL automorphism play a role?
For a directed path (v_0,...,v_m), the reversal is (v_m,...,v_0).
In a self-complementary tournament, reversal + complement gives a symmetry.
"""

import numpy as np
from itertools import permutations

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def enumerate_directed_paths(A, n, p):
    if p == 0: return [(v,) for v in range(n)]
    paths = []
    def dfs(path, depth):
        if depth == p:
            paths.append(tuple(path))
            return
        last = path[-1]
        visited = set(path)
        for v in range(n):
            if v not in visited and A[last][v]:
                path.append(v)
                dfs(path, depth+1)
                path.pop()
    for s in range(n):
        dfs([s], 0)
    return paths

# ============================================================
print("="*70)
print("PALEY SELF-COMPLEMENTARITY AND PATH REVERSAL")
print("="*70)

for p in [3, 5, 7]:
    if p % 4 != 3:
        print(f"\n  P_{p}: p ≡ {p%4} mod 4, NOT a Paley tournament")
        continue
    QR = {a*a % p for a in range(1, p)}
    NQR = set(range(1,p)) - QR
    A = circulant_tournament(p, QR)

    # Find a non-residue
    g = min(NQR)

    # The map x ↦ g*x mod p should be an anti-isomorphism
    # i.e., A_Paley[g*i mod p][g*j mod p] = 1 - A_Paley[i][j] for i≠j
    perm = [(g * i) % p for i in range(p)]
    anti_ok = True
    for i in range(p):
        for j in range(p):
            if i != j:
                if A[perm[i]][perm[j]] != 1 - A[i][j]:
                    anti_ok = False
                    break
    print(f"\n  P_{p}: QR={sorted(QR)}, NQR={sorted(NQR)}")
    print(f"    Non-residue g={g}, map x↦{g}x: anti-isomorphism? {anti_ok}")

    # For paths: the anti-automorphism maps directed m-paths of T
    # to directed m-paths of T^op.
    # Path reversal maps directed m-paths of T^op BACK to directed m-paths of T.
    # Combined: x ↦ g*x + reversal maps A_m(T) → A_m(T) but with some sign.

    # Count: how many directed m-paths are mapped to themselves?
    print(f"\n    Directed path counts:")
    for m in range(p):
        paths = enumerate_directed_paths(A, p, m)
        # Apply the anti-automorphism + reversal
        # This maps (v_0,...,v_m) with all v_i→v_{i+1} in T
        # to (g*v_m,...,g*v_0) with all g*v_{i+1}→g*v_i in T^op
        # and reversal gives (g*v_0,...,g*v_m) ... no wait.
        #
        # Anti-auto: (v_0,...,v_m) → (g*v_0,...,g*v_m) in T^op
        # Reversal: (g*v_0,...,g*v_m) → (g*v_m,...,g*v_0) in T
        # Combined: φ(v_0,...,v_m) = (g*v_m, g*v_{m-1}, ..., g*v_0) in T

        path_set = set(paths)
        fixed = 0
        for path in paths:
            img = tuple((g * path[m-i]) % p for i in range(m+1))
            if img == path:
                fixed += 1

        print(f"      m={m}: |A_m|={len(paths)}, fixed by φ={fixed}")

# ============================================================
print(f"\n{'='*70}")
print("PALEY P_7: DETAILED DUALITY ANALYSIS")
print("="*70)

p = 7
QR = {1,2,4}
NQR = {3,5,6}
A = circulant_tournament(p, QR)
g = 3  # smallest non-residue

print(f"  Anti-automorphism: x ↦ {g}x mod {p}")
print(f"  φ(v_0,...,v_m) = ({g}*v_m, ..., {g}*v_0) mod {p}")

# Check if φ induces a bijection A_m → A_{p-1-m}
# This would be the analogue of Poincaré duality
print(f"\n  Testing φ: A_m → A_{{p-1-m}} bijection:")

for m in range(p):
    paths_m = enumerate_directed_paths(A, p, m)
    paths_comp = enumerate_directed_paths(A, p, p-1-m)

    # φ: path ↦ complement of vertices used
    # For an m-path using m+1 vertices, the complement has p-m-1 vertices
    # = (p-1-m) + 1 - 1 = p-m-1 vertices in a (p-2-m)-path
    # This gives A_m → something of degree p-2-m, not p-1-m

    # Alternative: look at the relationship between Ω_m and Ω_{p-2-m}
    # P_7 Ω/7 = [1,3,6,9,9,6,3]
    # Indexing: m=0,1,2,3,4,5,6
    # p-2-m:   m=5,4,3,2,1,0,-1
    # This doesn't quite work either

    # Actually, the palindrome is Ω_m = Ω_{p-2-m}
    # m=0 ↔ m=5: Ω/7=[1] ↔ [6] — NOT equal (unless we consider something else)
    # Wait: [1,3,6,9,9,6,3] reversed = [3,6,9,9,6,3,1]
    # This IS the same for indices 1..5 but not 0 and 6

    pass

# Let me just verify the palindrome more carefully
print(f"\n  P_7 Ω/7 = [1, 3, 6, 9, 9, 6, 3]")
print(f"  Reversed = [3, 6, 9, 9, 6, 3, 1]")
print(f"  Ω_{0}/7 = 1, Ω_{6}/7 = 3  — NOT equal")
print(f"  Ω_{1}/7 = 3, Ω_{5}/7 = 6  — NOT equal")
print(f"  The palindrome claim was WRONG for P_7!")

# Wait, let me recheck P_5
print(f"\n  P_5 Ω/5 = [1, 2, 2, 2, 1]")
print(f"  This IS palindromic: Ω_0=Ω_4, Ω_1=Ω_3")

# P_7 nonzero Ω are [7, 21, 42, 63, 63, 42, 21]
# Not palindromic: 7 ≠ 21

# But earlier I said Ω_m/7 = [1,3,6,9,9,6,3]
# [1,3,6,9] and [9,6,3] — the sequence 9,6,3 mirrors 3,6,9
# So it's palindromic starting from index 1!

# Actually [1,3,6,9,9,6,3] reversed is [3,6,9,9,6,3,1]
# NOT the same. So NOT palindromic.

# Let me check what "near-palindromic" means
omega_7 = [1,3,6,9,9,6,3]
print(f"\n  Checking sub-palindromes:")
print(f"    [3,6,9,9,6,3] palindromic? {[3,6,9,9,6,3] == [3,6,9,9,6,3][::-1]}")
print(f"    The TAIL Ω_1..Ω_6 = [3,6,9,9,6,3] IS palindromic!")
print(f"    Only Ω_0 breaks the palindrome.")

# This makes sense: Ω_0 = p (vertices) is always special
# And Ω_{p-1} = |{Hamiltonian paths}| which is much smaller than Ω_0

# KEY: the sequence Ω_1, Ω_2, ..., Ω_{p-2} is palindromic for Paley!
# P_5: Ω_1..Ω_3 = [10, 10, 10] — trivially palindromic (constant)
# P_7: Ω_1..Ω_5 = [21, 42, 63, 63, 42] — palindromic!

omega_p5_mid = [10, 10, 10]
omega_p7_mid = [21, 42, 63, 63, 42]
print(f"\n  P_5: Ω_1..Ω_3 = {omega_p5_mid}, palindromic? {omega_p5_mid == omega_p5_mid[::-1]}")
print(f"  P_7: Ω_1..Ω_5 = {omega_p7_mid}, palindromic? {omega_p7_mid == omega_p7_mid[::-1]}")

# ============================================================
print(f"\n{'='*70}")
print("PATH COMPLEMENT MAP")
print("="*70)

# For a tournament T on [n], define the complement map:
# For a directed m-path σ = (v_0,...,v_m), the complement is
# the set of unused vertices V \ {v_0,...,v_m}.
# If all n vertices are used (m=n-1), the complement is empty.

# For GLMY on Paley P_p: we need a map Ω_m → Ω_{p-2-m}
# This requires a canonical way to turn the complement set into a path

# In the CLIQUE complex (complete graph), complement gives Poincaré duality.
# In the PATH complex, we need the complement vertices to form a path.

# For Paley, every subset of vertices induces a sub-tournament,
# and by self-complementarity, the complement sub-tournament is isomorphic.
# But do the complement vertices form a DIRECTED path?

# Let's check: for each Hamiltonian path (m=p-1), its complement is empty.
# For each (p-2)-path using p-1 vertices, the complement is 1 vertex.
# There are p * |A_{p-2} from 0| paths, and complement gives p vertices.

# More interesting: for m and p-2-m, do we get a bijection?
# Check at P_7: Ω_2 = 42, Ω_4 = 63. Not equal, so NOT a bijection.
# But Ω_1 = 21, Ω_5 = 42. Also not equal.

# Hmm. The complement map doesn't directly give a duality.
# Let me think about what DOES give the partial palindrome.

# ============================================================
print(f"\n{'='*70}")
print("RANK SEQUENCE ANALYSIS")
print("="*70)

# P_7 detailed:
# m:      0    1    2    3    4    5    6
# Ω_m:    7   21   42   63   63   42   21
# rk(∂):   -    6   15   27   36   21   21
# ker:     7   15   27   36   27   21   21
# im(+1):  6   15   27   36   21   21    0
# β:       1    0    0    0    6    0    0

# The ranks: rk(∂_m) = [_, 6, 15, 27, 36, 21, 21]
# rk/3 = [_, 2, 5, 9, 12, 7, 7]

# Note: rk(∂_5) = rk(∂_6) = 21. These are equal!
# And rk(∂_1) = 6 = Ω_0 - 1 = p - 1

# The β_4 = 6 = ker(∂_4) - im(∂_5) = (63-36) - 21 = 27 - 21 = 6

# So the key insight is: rk(∂_4) = 36 but rk(∂_5) = 21
# The jump from 36 to 21 creates the gap that gives β_4 = 6

# Is rk(∂_5)/rk(∂_4) related to something about Paley?
print(f"  P_7 rank ratios:")
rks = [None, 6, 15, 27, 36, 21, 21]
for m in range(2, 7):
    print(f"    rk(∂_{m})/rk(∂_{m-1}) = {rks[m]}/{rks[m-1]} = {rks[m]/rks[m-1]:.3f}")

# ============================================================
print(f"\n{'='*70}")
print("PALEY vs NON-PALEY: RANK COMPARISON")
print("="*70)

def compute_omega_basis(allowed_p, allowed_pm1, p_deg):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0,0))
    if p_deg <= 1: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for path in allowed_p:
        for i in range(p_deg+1):
            face = path[:i] + path[i+1:]
            if face not in allowed_pm1_set and face not in non_allowed:
                non_allowed[face] = na_count
                na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    J = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for i in range(p_deg+1):
            face = path[:i] + path[i+1:]
            if face in non_allowed:
                J[non_allowed[face], j] += (-1)**i
    U, s, Vh = np.linalg.svd(J, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    if rank == dim_Ap: return np.zeros((dim_Ap, 0))
    return Vh[rank:].T

def full_ranks(A):
    n = A.shape[0]
    allowed = {}
    omega = {}
    for m in range(n+1):
        allowed[m] = enumerate_directed_paths(A, n, m)
        omega[m] = compute_omega_basis(allowed[m], allowed.get(m-1, []), m)
    ranks = {}
    for m in range(1, n):
        dom = omega[m].shape[1]
        cod = omega[m-1].shape[1]
        if dom == 0 or cod == 0:
            ranks[m] = 0
            continue
        idx = {path: i for i, path in enumerate(allowed[m-1])}
        B = np.zeros((len(allowed[m-1]), len(allowed[m])))
        for j, path in enumerate(allowed[m]):
            for i in range(m+1):
                face = path[:i] + path[i+1:]
                if face in idx:
                    B[idx[face], j] += (-1)**i
        B_omega = omega[m-1].T @ B @ omega[m]
        ranks[m] = np.linalg.matrix_rank(B_omega, tol=1e-8)
    return {m: omega[m].shape[1] for m in range(n)}, ranks

for name, n_val, S in [("P_7 Paley", 7, {1,2,4}), ("I_7 Interval", 7, {1,2,3})]:
    A = circulant_tournament(n_val, S)
    om, rk = full_ranks(A)
    print(f"\n  {name}:")
    print(f"    Ω: {[om[m] for m in range(n_val)]}")
    print(f"    rk: {[rk.get(m, 0) for m in range(1, n_val)]}")

    # Compute β
    betti = []
    for m in range(n_val):
        b = om[m] - rk.get(m, 0) - rk.get(m+1, 0)
        betti.append(b)
    print(f"    β: {betti}")

print("\nDONE.")
