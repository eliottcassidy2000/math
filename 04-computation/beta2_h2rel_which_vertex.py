#!/usr/bin/env python3
"""
beta2_h2rel_which_vertex.py — Which vertices give h₂_rel=1 vs h₂_rel=0?

Understanding the pattern to prove HYP-258 algebraically.

Key questions:
1. What property of v determines h₂_rel?
2. Is it related to d+(v), β₁(T\v), or some other invariant?
3. Is there a simple criterion for h₂_rel = 0?

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def compute_h2_rel_details(A, n, v):
    """Compute h₂(T,T\\v) with full dimension details."""
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    ap0_T = enumerate_allowed_paths(A, n, 0)
    ap1_T = enumerate_allowed_paths(A, n, 1)
    ap2_T = enumerate_allowed_paths(A, n, 2)
    ap3_T = enumerate_allowed_paths(A, n, 3)

    ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
    ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
    ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
    ap3_sub = enumerate_allowed_paths(A_sub, n1, 3)

    om1_T = compute_omega_basis(A, n, 1, ap1_T, ap0_T)
    om2_T = compute_omega_basis(A, n, 2, ap2_T, ap1_T) if ap2_T else np.zeros((0,0))
    om3_T = compute_omega_basis(A, n, 3, ap3_T, ap2_T) if ap3_T else np.zeros((0,0))
    om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
    om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

    d2_T = dim_om(om2_T)
    d1_T = dim_om(om1_T)
    d3_T = dim_om(om3_T)
    d1_sub = dim_om(om1_sub)
    d2_sub = dim_om(om2_sub)

    # β₁(T\v)
    bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
    rk_d1 = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)
    z1_sub = d1_sub - rk_d1
    if d2_sub > 0:
        bd2_om = np.linalg.lstsq(om1_sub, build_full_boundary_matrix(ap2_sub, ap1_sub) @ om2_sub, rcond=None)[0]
        b1_sub = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    else:
        b1_sub = 0
    beta1_sub = z1_sub - b1_sub

    if d2_T == 0:
        return {'h2_rel': 0, 'beta1_sub': beta1_sub, 'd_rel2': 0, 'z2_rel': 0, 'b2_rel': 0}

    ap2_T_list = [tuple(p) for p in ap2_T]
    ap1_T_list = [tuple(p) for p in ap1_T]

    if ap2_sub and d2_sub > 0:
        embed = np.zeros((len(ap2_T_list), d2_sub))
        for j in range(d2_sub):
            for k, path_sub in enumerate(ap2_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap2_T_list:
                    embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
        phi = np.linalg.lstsq(om2_T, embed, rcond=None)[0]
    else:
        phi = np.zeros((d2_T, 0))

    rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
    if rk_phi > 0:
        U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
        Q = U_phi[:, rk_phi:]
    else:
        Q = np.eye(d2_T)

    d_rel = Q.shape[1]
    if d_rel == 0:
        return {'h2_rel': 0, 'beta1_sub': beta1_sub, 'd_rel2': 0, 'z2_rel': 0, 'b2_rel': 0}

    coords2_T = np.linalg.lstsq(om1_T, build_full_boundary_matrix(ap2_T, ap1_T) @ om2_T, rcond=None)[0]

    if d1_sub > 0:
        embed1 = np.zeros((len(ap1_T_list), d1_sub))
        for j in range(d1_sub):
            for k, path_sub in enumerate(ap1_sub):
                path_T = tuple(remap[x] for x in path_sub)
                if path_T in ap1_T_list:
                    embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
        psi = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]
    else:
        psi = np.zeros((d1_T, 0))

    rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
    if rk_psi > 0:
        U_psi, _, _ = np.linalg.svd(psi, full_matrices=True)
        R = U_psi[:, rk_psi:]
    else:
        R = np.eye(d1_T)

    coords2_rel = R.T @ coords2_T @ Q
    rk_d2_rel = np.linalg.matrix_rank(coords2_rel, tol=1e-8)
    z2_rel = d_rel - rk_d2_rel

    if z2_rel == 0:
        return {'h2_rel': 0, 'beta1_sub': beta1_sub, 'd_rel2': d_rel, 'z2_rel': 0, 'b2_rel': 0}

    if d3_T > 0:
        coords3_T = np.linalg.lstsq(om2_T, build_full_boundary_matrix(ap3_T, ap2_T) @ om3_T, rcond=None)[0]
        om3_sub = compute_omega_basis(A_sub, n1, 3, ap3_sub, ap2_sub) if ap3_sub else np.zeros((0,0))
        d3_sub = dim_om(om3_sub)
        d3_proj = Q.T @ coords3_T
        if d3_sub > 0:
            ap3_T_list = [tuple(p) for p in ap3_T]
            embed3 = np.zeros((len(ap3_T_list), d3_sub))
            for j in range(d3_sub):
                for k, path_sub in enumerate(ap3_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap3_T_list:
                        embed3[ap3_T_list.index(path_T), j] = om3_sub[k, j]
            chi = np.linalg.lstsq(om3_T, embed3, rcond=None)[0]
            rk_chi = np.linalg.matrix_rank(chi, tol=1e-8)
            if rk_chi > 0:
                U_chi, _, _ = np.linalg.svd(chi, full_matrices=True)
                d3_proj_rel = d3_proj @ U_chi[:, rk_chi:]
            else:
                d3_proj_rel = d3_proj
        else:
            d3_proj_rel = d3_proj
    else:
        d3_proj_rel = np.zeros((d_rel, 0))

    rk_d3_rel = np.linalg.matrix_rank(d3_proj_rel, tol=1e-8)
    b2_rel = rk_d3_rel
    h2_rel = z2_rel - b2_rel
    return {'h2_rel': h2_rel, 'beta1_sub': beta1_sub, 'd_rel2': d_rel, 'z2_rel': z2_rel, 'b2_rel': b2_rel}


print("=" * 70)
print("WHICH VERTEX GIVES h₂_rel = 0?")
print("=" * 70)

n = 5
m = n*(n-1)//2

# For each tournament, record which interior v give h2_rel=0 vs h2_rel=1
h2_by_dplus = defaultdict(lambda: Counter())  # d+ -> Counter of h2_rel values
h2_by_beta1 = defaultdict(lambda: Counter())
tours_pattern = Counter()  # pattern of which d+ give 0 vs 1

for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    interior_v = [v for v in range(n) if 1 <= scores[v] <= n-2]
    pattern = []

    for v in interior_v:
        d = compute_h2_rel_details(A, n, v)
        dv = scores[v]
        h2_by_dplus[dv][d['h2_rel']] += 1
        h2_by_beta1[d['beta1_sub']][d['h2_rel']] += 1
        pattern.append((dv, d['h2_rel'], d['beta1_sub']))

    # Sort pattern by d+
    pattern.sort()
    pattern_key = tuple((dv, h2r) for dv, h2r, _ in pattern)
    tours_pattern[pattern_key] += 1

print(f"\nn=5: h₂_rel distribution by d⁺(v):")
for dv in sorted(h2_by_dplus):
    print(f"  d⁺={dv}: {dict(h2_by_dplus[dv])}")

print(f"\nh₂_rel distribution by β₁(T\\v):")
for b1 in sorted(h2_by_beta1):
    print(f"  β₁={b1}: {dict(h2_by_beta1[b1])}")

print(f"\nTournament patterns (d⁺, h₂_rel):")
for pattern, cnt in sorted(tours_pattern.items(), key=lambda x: -x[1])[:20]:
    print(f"  {pattern}: {cnt}")

# KEY QUESTION: Is h₂_rel=0 ↔ β₁(T\v)=0?
print(f"\n{'='*70}")
print("KEY: Is h₂_rel correlated with β₁(T\\v)?")
print("=" * 70)

# Cross-tabulate h2_rel vs beta1_sub
cross = defaultdict(int)
for bits in range(1 << m):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    for v in range(n):
        dv = scores[v]
        if dv == 0 or dv == n-1:
            continue

        d = compute_h2_rel_details(A, n, v)
        cross[(d['h2_rel'], d['beta1_sub'])] += 1

print("\n  (h₂_rel, β₁(T\\v)) cross-tabulation:")
for (h2r, b1), cnt in sorted(cross.items()):
    print(f"    h₂_rel={h2r}, β₁={b1}: {cnt}")

# Is h2_rel=0 iff beta1_sub=0?
h2_0_b1_0 = cross.get((0, 0), 0)
h2_0_b1_pos = sum(cnt for (h2r, b1), cnt in cross.items() if h2r == 0 and b1 > 0)
h2_pos_b1_0 = sum(cnt for (h2r, b1), cnt in cross.items() if h2r > 0 and b1 == 0)
h2_pos_b1_pos = sum(cnt for (h2r, b1), cnt in cross.items() if h2r > 0 and b1 > 0)

print(f"\n  h₂_rel=0, β₁=0: {h2_0_b1_0}")
print(f"  h₂_rel=0, β₁>0: {h2_0_b1_pos}")
print(f"  h₂_rel>0, β₁=0: {h2_pos_b1_0}")
print(f"  h₂_rel>0, β₁>0: {h2_pos_b1_pos}")

if h2_pos_b1_0 == 0 and h2_0_b1_pos > 0:
    print(f"\n  h₂_rel>0 IMPLIES β₁>0 (necessary but not sufficient)")
    print(f"  Equivalently: β₁=0 IMPLIES h₂_rel=0")
    print(f"  So: if T\\v is acyclic (β₁=0), then h₂_rel=0 automatically!")
elif h2_0_b1_pos == 0:
    print(f"\n  h₂_rel=0 IFF β₁=0 (PERFECT CORRELATION)")
else:
    print(f"\n  No perfect correlation between h₂_rel and β₁")


# Does every tournament have some v with β₁(T\v) = 0?
print(f"\n{'='*70}")
print("Does every T have interior v with β₁(T\\v) = 0?")
print("=" * 70)

for n_test in [5, 6]:
    m_test = n_test*(n_test-1)//2
    total_test = 1 << m_test

    has_b1_zero = 0
    no_b1_zero = 0

    t0 = time.time()
    for bits in range(total_test):
        A = build_adj(n_test, bits)
        scores = [sum(A[i][j] for j in range(n_test) if j!=i) for i in range(n_test)]
        interior_v = [v for v in range(n_test) if 1 <= scores[v] <= n_test-2]

        found = False
        for v in interior_v:
            others = [i for i in range(n_test) if i != v]
            n1 = n_test - 1
            A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

            ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
            ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
            ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)

            om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
            om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

            d1_sub = dim_om(om1_sub)
            d2_sub = dim_om(om2_sub)

            bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
            rk_d1 = np.linalg.matrix_rank(bd1_sub @ om1_sub, tol=1e-8)
            z1 = d1_sub - rk_d1

            if d2_sub > 0:
                bd2_om = np.linalg.lstsq(om1_sub, build_full_boundary_matrix(ap2_sub, ap1_sub) @ om2_sub, rcond=None)[0]
                b1 = np.linalg.matrix_rank(bd2_om, tol=1e-8)
            else:
                b1 = 0

            beta1 = z1 - b1
            if beta1 == 0:
                found = True
                break

        if found:
            has_b1_zero += 1
        else:
            no_b1_zero += 1

        if (bits + 1) % 5000 == 0:
            print(f"  n={n_test}: {bits+1}/{total_test} ({time.time()-t0:.0f}s)")

    print(f"\nn={n_test}: {total_test} tournaments in {time.time()-t0:.0f}s")
    print(f"  Has interior v with β₁(T\\v)=0: {has_b1_zero}")
    print(f"  ALL interior v have β₁(T\\v)>0: {no_b1_zero}")

    if no_b1_zero == 0:
        print(f"  ✓ EVERY tournament has interior v with β₁(T\\v)=0")
    else:
        print(f"  ✗ {no_b1_zero} counterexamples!")


print("\nDone.")
