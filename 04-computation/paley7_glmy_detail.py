#!/usr/bin/env python3
"""
paley7_glmy_detail.py — opus-2026-03-13-S71

Deep dive into Paley P_7 GLMY path homology.
β = [1,0,0,0,6,0,0] — why β_4 = 6?

Also check: is β_4 = p-1 a pattern for Paley tournaments?
"""

import numpy as np

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

def compute_omega_basis(allowed_p, allowed_pm1, p):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0,0)), 0
    if p <= 1: return np.eye(dim_Ap), 0
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for path in allowed_p:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in allowed_pm1_set and face not in non_allowed:
                non_allowed[face] = na_count
                na_count += 1
    if na_count == 0: return np.eye(dim_Ap), 0
    J = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face in non_allowed:
                J[non_allowed[face], j] += (-1)**i
    U, s, Vh = np.linalg.svd(J, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    if rank == dim_Ap: return np.zeros((dim_Ap, 0)), na_count
    return Vh[rank:].T, na_count

def full_glmy(A):
    """Full GLMY with all intermediate data."""
    n = A.shape[0]
    max_dim = n - 1
    allowed = {}
    omega = {}
    junk_counts = {}
    for p in range(max_dim+2):
        allowed[p] = enumerate_directed_paths(A, n, p)
        ob, jc = compute_omega_basis(allowed[p], allowed.get(p-1, []), p)
        omega[p] = ob
        junk_counts[p] = jc

    ranks = {}
    kernels = {}
    images = {}
    for p in range(1, max_dim+1):
        dim_om = omega[p].shape[1]
        dim_om_prev = omega[p-1].shape[1]
        if dim_om == 0 or dim_om_prev == 0:
            ranks[p] = 0
            kernels[p] = dim_om
            images[p] = 0
            continue
        idx = {path: i for i, path in enumerate(allowed[p-1])}
        B = np.zeros((len(allowed[p-1]), len(allowed[p])))
        for j, path in enumerate(allowed[p]):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                if face in idx:
                    B[idx[face], j] += (-1)**i
        B_omega = omega[p-1].T @ B @ omega[p]
        rk = np.linalg.matrix_rank(B_omega, tol=1e-8)
        ranks[p] = rk
        kernels[p] = dim_om - rk
        images[p] = rk

    betti = []
    for p in range(max_dim+1):
        dim_p = omega[p].shape[1]
        rk_p = ranks.get(p, 0)
        rk_p1 = ranks.get(p+1, 0)
        betti.append(dim_p - rk_p - rk_p1)

    return {
        'betti': betti,
        'allowed': {p: len(allowed[p]) for p in range(max_dim+1)},
        'omega': {p: omega[p].shape[1] for p in range(max_dim+1)},
        'junk': junk_counts,
        'ranks': ranks,
        'kernels': kernels,
        'images': images,
    }

# ============================================================
print("="*70)
print("PALEY P_7 GLMY DETAILED ANALYSIS")
print("="*70)

QR7 = {1,2,4}
A = circulant_tournament(7, QR7)
data = full_glmy(A)

print(f"\n  Betti:   {data['betti']}")
print(f"  |A_m|:   {[data['allowed'][p] for p in sorted(data['allowed'])]}")
print(f"  dim(Ω):  {[data['omega'][p] for p in sorted(data['omega'])]}")
print(f"  #junk:   {[data['junk'][p] for p in sorted(data['junk'])]}")
print(f"  rk(∂):   {[data['ranks'].get(p,0) for p in range(1,7)]}")
print(f"  ker(∂):  {[data['kernels'].get(p,'-') for p in range(1,7)]}")
print(f"  im(∂):   {[data['images'].get(p,'-') for p in range(1,7)]}")

print(f"\n  Chain: Ω_6→Ω_5→Ω_4→Ω_3→Ω_2→Ω_1→Ω_0")
om = [data['omega'][p] for p in range(7)]
rk = [data['ranks'].get(p,0) for p in range(1,8)]
for p in range(7):
    ker = om[p] - rk[p] if p > 0 else om[p]
    im_into = rk[p+1] if p+1 <= 6 else 0
    beta = data['betti'][p]
    print(f"  Ω_{p}={om[p]:3d}  rk(∂_{p})={rk[p]:3d}  ker(∂_{p})={ker:3d}  im(∂_{p+1})={im_into:3d}  β_{p}={beta}")

# ============================================================
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC AND PALEY PATTERN")
print("="*70)

# Check Paley at n=3, 5, 7
for p, QR in [(3, {1}), (5, {1,2}), (7, {1,2,4})]:
    A = circulant_tournament(p, QR)
    d = full_glmy(A)
    chi = sum((-1)**m * d['betti'][m] for m in range(len(d['betti'])))
    print(f"  Paley P_{p}: β={d['betti']}, chi={chi}")
    if chi != 0:
        print(f"    chi/p = {chi/p}")

# Check if there's a pattern: β_{p-3}(P_p) = p-1
print(f"\n  Pattern check: β_{{p-3}}(P_p) = p-1?")
for p, QR in [(3, {1}), (5, {1,2}), (7, {1,2,4})]:
    A = circulant_tournament(p, QR)
    d = full_glmy(A)
    idx = p - 3
    if 0 <= idx < len(d['betti']):
        print(f"    P_{p}: β_{idx} = {d['betti'][idx]}, p-1 = {p-1}")

# ============================================================
print(f"\n{'='*70}")
print("C7 INTERVAL vs PALEY: WHY DIFFERENT β?")
print("="*70)

for name, S in [("Paley {1,2,4}", {1,2,4}), ("Interval {1,2,3}", {1,2,3}),
                ("Balanced {1,3,5}", {1,3,5})]:
    A = circulant_tournament(7, S)
    d = full_glmy(A)
    print(f"\n  {name}:")
    for p in range(7):
        print(f"    m={p}: |A|={d['allowed'][p]:4d}, Ω={d['omega'][p]:3d}, "
              f"junk_faces={d['junk'][p]:4d}, gap={d['allowed'][p]-d['omega'][p]:4d}")

# ============================================================
print(f"\n{'='*70}")
print("GLMY Ω/n SEQUENCES FOR CIRCULANTS")
print("="*70)

print("  Hypothesis: Ω_m/n for circulant on Z_n always integer")
print()

cases = [
    ("C3", 3, {1}),
    ("C5", 5, {1,2}),
    ("P7", 7, {1,2,4}),
    ("I7", 7, {1,2,3}),
    ("B7", 7, {1,3,5}),
]

for name, n, S in cases:
    A = circulant_tournament(n, S)
    d = full_glmy(A)
    om = [d['omega'][p] for p in range(n)]
    om_n = [o//n for o in om]
    print(f"  {name}: Ω/n = {om_n}")

# ============================================================
# Ω/n for Paley is [1,3,6,9,9,6,3] — check OEIS or binomial pattern
print(f"\n  Paley P_7 Ω/n: [1,3,6,9,9,6,3]")
print(f"  Binomial C(6,m): [1,6,15,20,15,6,1]")
print(f"  C(5,m):          [1,5,10,10,5,1]")
print(f"  Note: sum(Ω/n) = {sum([1,3,6,9,9,6,3])}")
print(f"  sum(C(6,m))     = {sum([1,6,15,20,15,6,1])}")

# Cumulative sums / partial sums?
seq = [1,3,6,9,9,6,3]
diffs = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
print(f"  Differences: {diffs}")
print(f"  Second diff: {[diffs[i+1]-diffs[i] for i in range(len(diffs)-1)]}")

# Check: is this related to (n-1 choose m) somehow?
# At n=7: q=(7+1)/4=2, Q-1=1
# The Ω/n for Paley are triangular-like: 1, 1+2, 1+2+3, ...
# 1,3,6 = C(m+1,2) for m=0,1,2
# Then 9,9,6,3 ...
from math import comb
print(f"  C(m+1,2): {[comb(m+1,2) for m in range(7)]}")
print(f"  C(m+2,3): {[comb(m+2,3) for m in range(7)]}")

# ============================================================
print(f"\n{'='*70}")
print("DIRECTED PATH COUNTS / OMEGA STRUCTURE")
print("="*70)

# Key insight: for Paley P_7, GLMY Ω_4=63 >> TRH regular_4=21
# The GLMY uses ALL directed paths, not just regular ones
# The Ω subspace filters via junk cancellation, giving very different answers

print("  For P_7:")
print("    Regular paths use: v_i→v_{i+1} AND v_{i-1}→v_{i+1}")
print("    Directed paths use: v_i→v_{i+1} only")
print("    Regular ⊂ Directed always")
print()

A = circulant_tournament(7, {1,2,4})
for m in range(7):
    dp = enumerate_directed_paths(A, 7, m)
    # Count how many directed paths are also regular
    regular_count = 0
    for path in dp:
        is_reg = True
        for i in range(1, len(path)-1):
            if not A[path[i-1]][path[i+1]]:
                is_reg = False
                break
        regular_count += 1 if is_reg else 0
    print(f"  m={m}: directed={len(dp):4d}, regular_in_directed={regular_count:4d}")

print("\nDONE.")
