#!/usr/bin/env python3
"""
paley_correct_betti.py — opus-2026-03-13-S71

Correct Paley tournament GLMY Betti numbers.
Paley tournament P_p is defined for p ≡ 3 mod 4:
  i→j iff j-i is a quadratic residue mod p.

For p ≡ 1 mod 4, QR = -QR so this gives a symmetric (non-tournament) graph.

Valid Paley tournaments: p = 3, 7, 11, 19, 23, ...
"""

import numpy as np

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def verify_tournament(A):
    n = A.shape[0]
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] + A[j][i] != 1:
                return False
    return True

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

def compute_omega_basis(allowed_p, allowed_pm1, p):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0,0))
    if p <= 1: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for path in allowed_p:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in allowed_pm1_set and face not in non_allowed:
                non_allowed[face] = na_count
                na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    J = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face in non_allowed:
                J[non_allowed[face], j] += (-1)**i
    U, s, Vh = np.linalg.svd(J, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    if rank == dim_Ap: return np.zeros((dim_Ap, 0))
    return Vh[rank:].T

def glmy_betti(A, max_dim=None):
    n = A.shape[0]
    if max_dim is None:
        max_dim = n - 1
    allowed = {}
    omega = {}
    for p in range(max_dim+2):
        allowed[p] = enumerate_directed_paths(A, n, p)
        omega[p] = compute_omega_basis(allowed[p], allowed.get(p-1, []), p)
    ranks = {}
    for p in range(1, max_dim+1):
        dim_om = omega[p].shape[1]
        dim_om_prev = omega[p-1].shape[1]
        if dim_om == 0 or dim_om_prev == 0:
            ranks[p] = 0
            continue
        idx = {path: i for i, path in enumerate(allowed[p-1])}
        B = np.zeros((len(allowed[p-1]), len(allowed[p])))
        for j, path in enumerate(allowed[p]):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                if face in idx:
                    B[idx[face], j] += (-1)**i
        B_omega = omega[p-1].T @ B @ omega[p]
        ranks[p] = np.linalg.matrix_rank(B_omega, tol=1e-8)
    betti = []
    for p in range(max_dim+1):
        dim_p = omega[p].shape[1]
        rk_p = ranks.get(p, 0)
        rk_p1 = ranks.get(p+1, 0)
        betti.append(dim_p - rk_p - rk_p1)
    omega_dims = [omega[p].shape[1] for p in range(max_dim+1)]
    a_dims = [len(allowed[p]) for p in range(max_dim+1)]
    return betti, omega_dims, a_dims

# ============================================================
print("="*70)
print("PALEY TOURNAMENT VERIFICATION")
print("="*70)

for p in [3, 5, 7, 11, 13, 19, 23]:
    QR = {a*a % p for a in range(1, p)}
    neg_QR = {p - q for q in QR}
    is_tournament = (QR & neg_QR == set())
    print(f"  p={p:2d}: p mod 4 = {p%4}, QR={sorted(QR)}, -QR={sorted(neg_QR)}, "
          f"tournament? {is_tournament}, QR=-QR? {QR == neg_QR}")

# ============================================================
print(f"\n{'='*70}")
print("CORRECT PALEY TOURNAMENT GLMY BETTI")
print("="*70)

for p in [3, 7]:
    QR = {a*a % p for a in range(1, p)}
    A = circulant_tournament(p, QR)
    assert verify_tournament(A), f"P_{p} is not a tournament!"

    b, o, a = glmy_betti(A)
    chi = sum((-1)**m * b[m] for m in range(len(b)))

    print(f"\n  Paley P_{p}: QR={sorted(QR)}")
    print(f"    A_m:  {a}")
    print(f"    Ω_m:  {o}")
    print(f"    β_m:  {b}")
    print(f"    chi = {chi}")
    if all(o_m % p == 0 for o_m in o):
        print(f"    Ω/p:  {[x//p for x in o]}")

# ============================================================
print(f"\n{'='*70}")
print("ALL n=5 CIRCULANT TOURNAMENT TYPES")
print("="*70)

n = 5
print(f"\n  Valid orientation sets S for n={n} (S ∩ -S = ∅, S ∪ -S = {{1..{n-1}}}):")
from itertools import combinations
valid_S = []
for S_tuple in combinations(range(1, n), (n-1)//2):
    S = set(S_tuple)
    neg_S = {n - s for s in S}
    if S & neg_S == set() and S | neg_S == set(range(1, n)):
        valid_S.append(S)
        A = circulant_tournament(n, S)
        b, o, a = glmy_betti(A)
        chi = sum((-1)**m * b[m] for m in range(len(b)))
        print(f"    S={sorted(S)}: β={b}, Ω={o}, chi={chi}")

# ============================================================
print(f"\n{'='*70}")
print("ALL n=7 CIRCULANT TOURNAMENT TYPES")
print("="*70)

n = 7
valid_S = []
for S_tuple in combinations(range(1, n), (n-1)//2):
    S = set(S_tuple)
    neg_S = {n - s for s in S}
    if S & neg_S == set() and S | neg_S == set(range(1, n)):
        valid_S.append(S)

from itertools import permutations as perms

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in perms(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best: best = enc
    return best

# Group by isomorphism type
iso_classes = {}
for S in valid_S:
    A = circulant_tournament(n, S)
    cf = canon_form(A)
    if cf not in iso_classes:
        iso_classes[cf] = []
    iso_classes[cf].append(sorted(S))

print(f"  {len(valid_S)} valid orientation sets, {len(iso_classes)} isomorphism types")

for cf, Slist in sorted(iso_classes.items()):
    A = circulant_tournament(n, set(Slist[0]))
    b, o, a = glmy_betti(A)
    chi = sum((-1)**m * b[m] for m in range(len(b)))
    print(f"\n  Type (represented by S={Slist[0]}):")
    print(f"    All S in this class: {Slist}")
    print(f"    β = {b}, chi = {chi}")
    print(f"    Ω = {o}")
    if all(x % n == 0 for x in o):
        print(f"    Ω/n = {[x//n for x in o]}")

# ============================================================
print(f"\n{'='*70}")
print("PATTERN ANALYSIS: PALEY vs NON-PALEY CIRCULANTS")
print("="*70)

# Summary table
print(f"\n  {'Name':20s} {'β':40s} {'chi':5s} {'Ω/n':30s}")
print(f"  {'-'*20} {'-'*40} {'-'*5} {'-'*30}")

for name, n_val, S in [
    ("C3 (Paley)", 3, {1}),
    ("C5_{1,2}", 5, {1,2}),
    ("C7 Paley {1,2,4}", 7, {1,2,4}),
    ("C7 Interval {1,2,3}", 7, {1,2,3}),
]:
    A = circulant_tournament(n_val, S)
    b, o, a = glmy_betti(A)
    chi = sum((-1)**m * b[m] for m in range(len(b)))
    o_n = [x//n_val for x in o] if all(x % n_val == 0 for x in o) else "N/A"
    print(f"  {name:20s} {str(b):40s} {chi:5d} {str(o_n):30s}")

print("\nDONE.")
