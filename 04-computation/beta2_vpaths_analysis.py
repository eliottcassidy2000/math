#!/usr/bin/env python3
"""
beta2_vpaths_analysis.py — Detailed structure of relative Ω₂ for vertex v

The relative chain space Ω₂^rel = Ω₂(T)/Ω₂(T\\v) is spanned by
those Ω₂-elements that involve vertex v.

For a 2-path (a,b,c) to be allowed: a→b and b→c (edges exist).
It's in Ω₂ if ∂₁∂₂(a,b,c) = 0, i.e., (a,c) - (a,b) - (b,c) + ... = 0 in Ω₀.

A "v-path" is a path containing v. Three positions: v=a (start), v=b (mid), v=c (end).

Key question: what structural property of v determines dim(Ω₂^rel)?

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


print("=" * 70)
print("V-PATH STRUCTURE IN RELATIVE Ω₂")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For each interior vertex with h₂_rel=1, analyze the Ω₂ generator
# that represents the relative H₂ class

# First: understand allowed 2-paths through v by position
print(f"\n--- v at each position in allowed 2-paths ---")

position_data = defaultdict(list)

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]
    ap2 = enumerate_allowed_paths(A, n, 2)

    for v in range(n):
        if scores[v] == 0 or scores[v] == n - 1:
            continue

        v_start = sum(1 for p in ap2 if p[0] == v)  # v at position 0
        v_mid = sum(1 for p in ap2 if p[1] == v)    # v at position 1
        v_end = sum(1 for p in ap2 if p[2] == v)    # v at position 2

        position_data[(scores[v], v_start, v_mid, v_end)].append(bits)

print(f"(d+, #start, #mid, #end) → count of (tournament, vertex) pairs:")
for key in sorted(position_data, key=lambda k: (k[0], -len(position_data[k]))):
    print(f"  d+={key[0]}, pos=(start={key[1]}, mid={key[2]}, end={key[3]}): {len(position_data[key])} pairs")


# Now: for each (tournament, interior vertex), compute:
# 1. The Ω₂ basis vectors involving v
# 2. Whether they form non-transitive triples
# 3. h₂_rel value
print(f"\n--- Ω₂ basis analysis for v-elements ---")

# Simpler approach: for each tournament, look at ALL Ω₂ basis vectors,
# classify which ones involve v, and see how the relative dimension arises

n = 5
results = defaultdict(list)

for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap2:
        continue

    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    d2 = dim_om(om2)

    for v in range(n):
        if scores[v] == 0 or scores[v] == n-1:
            continue

        # Which allowed 2-paths involve v?
        v_indices = [i for i, p in enumerate(ap2) if v in p]
        non_v_indices = [i for i, p in enumerate(ap2) if v not in p]

        # The Ω₂ basis in coordinates of allowed 2-paths
        # om2[i,j] = coefficient of i-th allowed 2-path in j-th basis vector
        # "involves v" = has nonzero entry on v_indices rows

        # Relative Ω₂: quotient by sub-Ω₂ (non-v part)
        # In practice: compute rank of om2 restricted to non-v rows
        if non_v_indices:
            sub_block = om2[non_v_indices, :]
            rk_sub = np.linalg.matrix_rank(sub_block, tol=1e-8)
        else:
            rk_sub = 0

        d_rel = d2 - rk_sub  # dimension of relative Ω₂

        # Now compute actual h₂_rel
        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

        d2_sub = dim_om(om2_sub)

        # h₂_rel via the standard method
        remap = {i: others[i] for i in range(n1)}
        ap2_T_list = [tuple(p) for p in ap2]

        if ap2_sub and d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2, embed, rcond=None)[0]
            rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        else:
            rk_phi = 0

        # Check: d_rel should equal d2 - rk_phi
        d_rel_check = d2 - rk_phi

        results[(d2, d2_sub, rk_phi, d_rel_check, rk_sub, d_rel)].append((bits, v, scores[v]))

print(f"(d2_T, d2_sub, rk_incl, d_rel_incl, rk_sub_nonv, d_rel_nonv) → count:")
for key in sorted(results, key=lambda k: -len(results[k])):
    cnt = len(results[key])
    print(f"  {key} → {cnt}  {'✓' if key[3]==key[5] else '✗ MISMATCH'}")


# Part 3: For h₂_rel=1 cases, what is the actual relative cycle?
print(f"\n--- Part 3: Structure of relative H₂ generator ---")
print("For first 5 cases with h₂_rel=1:")

count = 0
for bits in range(total):
    if count >= 5:
        break
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]

    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    if not ap2:
        continue

    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))
    d2 = dim_om(om2)
    d1 = dim_om(om1)
    d3 = dim_om(om3)

    for v in range(n):
        if count >= 5:
            break
        if scores[v] == 0 or scores[v] == n-1:
            continue

        others = [i for i in range(n) if i != v]
        n1 = n - 1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        remap = {i: others[i] for i in range(n1)}
        ap2_T_list = [tuple(p) for p in ap2]
        ap1_T_list = [tuple(p) for p in ap1]

        ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)
        ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
        ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
        om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))
        om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
        d2_sub = dim_om(om2_sub)
        d1_sub = dim_om(om1_sub)

        if d2_sub > 0:
            embed = np.zeros((len(ap2_T_list), d2_sub))
            for j in range(d2_sub):
                for k, path_sub in enumerate(ap2_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap2_T_list:
                        embed[ap2_T_list.index(path_T), j] = om2_sub[k, j]
            phi = np.linalg.lstsq(om2, embed, rcond=None)[0]
        else:
            phi = np.zeros((d2, 0))

        rk_phi = np.linalg.matrix_rank(phi, tol=1e-8)
        if rk_phi > 0:
            U_phi, _, _ = np.linalg.svd(phi, full_matrices=True)
            Q = U_phi[:, rk_phi:]
        else:
            Q = np.eye(d2)

        d_rel = Q.shape[1]
        if d_rel == 0:
            continue

        # Compute relative ∂₂
        coords2_T = np.linalg.lstsq(om1, build_full_boundary_matrix(ap2, ap1) @ om2, rcond=None)[0]

        if d1_sub > 0:
            embed1 = np.zeros((len(ap1_T_list), d1_sub))
            for j in range(d1_sub):
                for k, path_sub in enumerate(ap1_sub):
                    path_T = tuple(remap[x] for x in path_sub)
                    if path_T in ap1_T_list:
                        embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
            psi = np.linalg.lstsq(om1, embed1, rcond=None)[0]
        else:
            psi = np.zeros((d1, 0))

        rk_psi = np.linalg.matrix_rank(psi, tol=1e-8)
        R = np.linalg.svd(psi, full_matrices=True)[0][:, rk_psi:] if rk_psi > 0 else np.eye(d1)

        coords2_rel = R.T @ coords2_T @ Q
        rk_d2_rel = np.linalg.matrix_rank(coords2_rel, tol=1e-8)
        z2_rel = d_rel - rk_d2_rel
        if z2_rel == 0:
            continue

        # This vertex has h₂_rel = 1. Find the relative cycle.
        # Kernel of coords2_rel
        U, S, Vt = np.linalg.svd(coords2_rel)
        # Null space of coords2_rel: last columns of Vt.T
        null_vec = Vt[-1, :]  # in Q-coordinates

        # Convert to Ω₂ coordinates
        omega_vec = Q @ null_vec
        # Convert to allowed-2-path coordinates
        path_vec = om2 @ omega_vec

        # Show which paths contribute
        print(f"\nTournament bits={bits}, v={v} (d+={scores[v]}), h₂_rel=1:")
        print(f"  Adjacency: {[''.join(str(A[i][j]) for j in range(n)) for i in range(n)]}")
        paths_with_coeff = []
        for i, p in enumerate(ap2):
            if abs(path_vec[i]) > 1e-8:
                # Is this path transitive? a→c exists?
                tt = "TT" if A[p[0]][p[2]] else "NT"
                pos = "start" if p[0]==v else ("mid" if p[1]==v else "end")
                paths_with_coeff.append((p, path_vec[i], tt, pos))

        for p, c, tt, pos in paths_with_coeff:
            print(f"  {c:+.3f} * {tuple(p)} [{tt}, v@{pos}]")

        # Check: what's the boundary of this element?
        bd_vec = build_full_boundary_matrix(ap2, ap1) @ path_vec
        # Express in terms of allowed 1-paths
        bd_norm = np.linalg.norm(bd_vec)
        print(f"  ∂₂ norm = {bd_norm:.6f}")
        if bd_norm > 1e-8:
            for i, p in enumerate(ap1):
                if abs(bd_vec[i]) > 1e-8:
                    v_in = "v-path" if v in p else "non-v"
                    print(f"    ∂: {bd_vec[i]:+.3f} * {tuple(p)} [{v_in}]")

        count += 1

print("\nDone.")
