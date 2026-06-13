#!/usr/bin/env python3
"""
beta2_h1_killing.py — H₁-killing mechanism: how vertex v fills cycles

Reformulation: h₂_rel(T,T\\v) = 1 iff v "fills" a cycle in T\\v.
More precisely: ∃ cycle class z ∈ H₁(T\\v) killed by inclusion into T,
AND this killing provides a nonzero element in ker(δ) of the LES.

The question: what does the "filling" look like combinatorially?

For each (T, v) with h₂_rel = 1:
1. Find the H₁ generator of T\\v (a 1-cycle)
2. Show it's in im(∂₂) in T (it becomes a boundary)
3. The ∂₂ preimage must involve v-paths

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
print("H₁-KILLING MECHANISM: WHAT VERTEX v DOES TO CYCLES")
print("=" * 70)

n = 5

# Study the first few tournaments with Σ h₂_rel = 3
for bits in range(1024):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    # Quick check: this tournament has the right score sequence for Σ=3
    # From earlier: Σ=3 tournaments have t₃=4 at n=5

    # Compute t₃
    t3 = 0
    cycles_list = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    t3 += 1
                    cycles_list.append((i,j,k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    t3 += 1
                    cycles_list.append((i,k,j))

    if t3 != 4:
        continue

    # Compute beta1(T)
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d1 = dim_om(om1)
    d2 = dim_om(om2)
    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk1 = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    z1 = d1 - rk1
    if d2 > 0:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2om = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        b1 = np.linalg.matrix_rank(bd2om, tol=1e-8)
    else:
        b1 = 0
    beta1_T = z1 - b1

    if beta1_T > 0:
        continue  # HYP-263: β₁>0 ⟹ Σ=0

    # This is a β₁=0, t₃=4 tournament — should have Σ=3
    print(f"\nTournament bits={bits}, t₃={t3}, β₁={beta1_T}")
    print(f"  Scores: {scores}")
    adj = [''.join(str(A[i][j]) for j in range(n)) for i in range(n)]
    print(f"  Adj: {adj}")

    # List all 4 3-cycles
    print(f"  3-cycles: {cycles_list}")

    # For each vertex v: which cycles survive in T\\v? What's β₁(T\\v)?
    for v in range(n):
        others = [i for i in range(n) if i != v]
        surviving_cycles = [c for c in cycles_list if v not in c]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

        # β₁(T\\v)
        ap0s = enumerate_allowed_paths(A_sub, n1, 0)
        ap1s = enumerate_allowed_paths(A_sub, n1, 1)
        ap2s = enumerate_allowed_paths(A_sub, n1, 2)
        om1s = compute_omega_basis(A_sub, n1, 1, ap1s, ap0s)
        om2s = compute_omega_basis(A_sub, n1, 2, ap2s, ap1s) if ap2s else np.zeros((0,0))
        d1s = dim_om(om1s)
        bd1s = build_full_boundary_matrix(ap1s, ap0s)
        rk1s = np.linalg.matrix_rank(bd1s @ om1s, tol=1e-8)
        z1s = d1s - rk1s
        if dim_om(om2s) > 0:
            bd2s = build_full_boundary_matrix(ap2s, ap1s)
            bd2oms = np.linalg.lstsq(om1s, bd2s @ om2s, rcond=None)[0]
            b1s = np.linalg.matrix_rank(bd2oms, tol=1e-8)
        else:
            b1s = 0
        beta1_sub = z1s - b1s

        # Count cycles through v
        cyc_v = sum(1 for c in cycles_list if v in c)

        print(f"  v={v} (d+={scores[v]}): cyc(v)={cyc_v}, surviving={len(surviving_cycles)}, β₁(T\\v)={beta1_sub}")

        if beta1_sub > 0:
            # Find the H₁ generator of T\\v
            # Z₁ basis in Ω₁(T\\v)
            U_z, S_z, Vt_z = np.linalg.svd(bd1s @ om1s, full_matrices=True)
            rk = sum(s > 1e-8 for s in S_z)
            z1_basis = Vt_z[rk:, :]  # Z₁ basis in Ω₁ coords

            # H₁ generator: Z₁ mod B₁
            if b1s > 0:
                b1_in_z1 = z1_basis @ bd2oms
                U_b, S_b, _ = np.linalg.svd(b1_in_z1, full_matrices=True)
                b_rk = sum(s > 1e-8 for s in S_b)
                h1_gen_z = U_b[:, b_rk:][:, 0]  # in Z₁ coords
            else:
                h1_gen_z = z1_basis[0]  # first Z₁ basis vector

            # Convert to A₁ coords
            h1_gen_sub = om1s @ (z1_basis.T @ h1_gen_z)  # in A₁(T\\v)

            # Show which edges are in the H₁ generator
            remap = {i: others[i] for i in range(n1)}
            gen_edges = []
            for i, p in enumerate(ap1s):
                if abs(h1_gen_sub[i]) > 1e-8:
                    gen_edges.append((remap[p[0]], remap[p[1]], h1_gen_sub[i]))
            print(f"    H₁ gen: {[(a,b,round(c,3)) for a,b,c in gen_edges]}")

            # Now: does this H₁ class become a boundary in T?
            # Embed the cycle into A₁(T)
            ap1_T_list = [tuple(p) for p in ap1]
            h1_in_T = np.zeros(len(ap1_T_list))
            for i, p in enumerate(ap1s):
                path_T = (remap[p[0]], remap[p[1]])
                if path_T in ap1_T_list:
                    h1_in_T[ap1_T_list.index(path_T)] = h1_gen_sub[i]

            # Express in Ω₁(T) coordinates
            h1_om_T = np.linalg.lstsq(om1, h1_in_T, rcond=None)[0]

            # Check: is this in B₁(T) = im(∂₂)?
            if d2 > 0:
                resid = h1_om_T - bd2om @ np.linalg.lstsq(bd2om, h1_om_T, rcond=None)[0]
                resid_norm = np.linalg.norm(resid)
                print(f"    In B₁(T)? residual = {resid_norm:.6f} ({'YES' if resid_norm < 1e-6 else 'NO'})")

                if resid_norm < 1e-6:
                    # Find the ∂₂ preimage (which Ω₂ elements fill this cycle)
                    preimage = np.linalg.lstsq(bd2om, h1_om_T, rcond=None)[0]
                    # This is in Ω₂(T) coordinates. Convert to A₂ paths.
                    preimage_A2 = om2 @ preimage
                    fill_paths = []
                    for i, p in enumerate(ap2):
                        if abs(preimage_A2[i]) > 1e-8:
                            has_v = v in p
                            tt = A[p[0]][p[2]]
                            fill_paths.append((tuple(p), round(preimage_A2[i], 3), has_v, 'TT' if tt else 'NT'))
                    print(f"    Filling 2-chain: {len(fill_paths)} paths")
                    for p, c, hv, tt in fill_paths[:10]:
                        print(f"      {c:+.3f} * {p} [{'v-path' if hv else 'non-v'}, {tt}]")
            else:
                print(f"    B₁(T) = 0, not filled")

    # Only show one tournament
    break

print("\nDone.")
