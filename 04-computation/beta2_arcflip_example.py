#!/usr/bin/env python3
"""
beta2_arcflip_example.py - Concrete examples of arc-flip invariance mechanism

Find specific arc flips where dOm2=0, db1=+1 (hardest non-trivial case):
- Same number of TT triples before and after
- beta_1 increases by 1 (a 1-cycle "appears")
- rk(d_3) also increases by 1 (a 2-cycle gets "filled")

Explicitly identify the new 1-cycle and the new 2-boundary.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def flip_arc(bits, i, j, n):
    idx = 0
    for a in range(n):
        for b in range(a+1, n):
            if a == i and b == j:
                return bits ^ (1 << idx)
            idx += 1
    return bits

def compute_chain_data(A, n):
    """Compute full chain complex data."""
    paths = {}
    omega = {}
    for p in range(5):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    # Boundary matrices
    bds = {}
    for p in [1, 2, 3, 4]:
        if len(paths.get(p, [])) > 0 and len(paths.get(p-1, [])) > 0:
            bds[p] = build_full_boundary_matrix(paths[p], paths[p-1])
        else:
            bds[p] = None

    # Boundary in Omega coordinates
    bd_om = {}
    for p in [2, 3]:
        dim_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_p > 0 and bds.get(p) is not None:
            raw = bds[p] @ omega[p]
            if p >= 2:
                bd_om[p], _, _, _ = np.linalg.lstsq(omega[p-1], raw, rcond=None)
            else:
                bd_om[p] = raw
        else:
            bd_om[p] = None

    # Ranks
    rks = {}
    for p in [1, 2, 3]:
        if bd_om.get(p) is not None:
            rks[p] = np.linalg.matrix_rank(bd_om[p], tol=1e-8)
        elif bds.get(p) is not None and omega.get(p) is not None:
            dim_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
            if dim_p > 0:
                raw = bds[p] @ omega[p]
                rks[p] = np.linalg.matrix_rank(raw, tol=1e-8)
            else:
                rks[p] = 0
        else:
            rks[p] = 0

    # Betti numbers
    dims = {p: (omega[p].shape[1] if omega[p].ndim == 2 else 0) for p in range(5)}
    betas = {}
    for p in range(4):
        z_p = dims[p] - rks.get(p, 0)
        b_p = z_p - rks.get(p+1, 0)
        betas[p] = b_p

    return paths, omega, bds, bd_om, dims, rks, betas


# Find a specific example at n=5
n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print("FINDING CONCRETE ARC-FLIP EXAMPLE (dOm2=0, db1=+1)")
print("=" * 70)

found = False
for bits in range(total):
    A = build_adj(n, bits)
    p1, o1, bd1, bdo1, d1, r1, b1 = compute_chain_data(A, n)

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            A2 = build_adj(n, bits2)
            p2, o2, bd2, bdo2, d2, r2, b2 = compute_chain_data(A2, n)

            dOm2 = d2[2] - d1[2]
            db1 = b2[1] - b1[1]
            drk3 = r2[3] - r1[3]

            if dOm2 == 0 and db1 == 1:
                if A[i][j] == 1:
                    u, v = i, j
                else:
                    u, v = j, i

                print(f"\nExample: bits={bits}, flip arc ({u}->{v}) to ({v}->{u})")
                print(f"  Scores T: {[sum(row) for row in A]}")
                print(f"  Scores T': {[sum(row) for row in A2]}")
                print(f"  beta(T): {[b1[p] for p in range(4)]}")
                print(f"  beta(T'): {[b2[p] for p in range(4)]}")
                print(f"  dim(Om2): {d1[2]} -> {d2[2]} (delta={dOm2})")
                print(f"  rk(d2): {r1[2]} -> {r2[2]} (delta={r2[2]-r1[2]})")
                print(f"  dim(Z2): {d1[2]-r1[2]} -> {d2[2]-r2[2]}")
                print(f"  rk(d3): {r1[3]} -> {r2[3]} (delta={drk3})")

                # List TT triples in T and T'
                tt_T = [(a,b,c) for (a,b,c) in p1[2]
                        if A[a][b] and A[b][c] and A[a][c]]
                tt_T2 = [(a,b,c) for (a,b,c) in p2[2]
                         if A2[a][b] and A2[b][c] and A2[a][c]]

                lost = [t for t in tt_T if t not in tt_T2]
                gained = [t for t in tt_T2 if t not in tt_T]
                stable = [t for t in tt_T if t in tt_T2]

                print(f"\n  TT triples T: {len(tt_T)}")
                print(f"  TT triples T': {len(tt_T2)}")
                print(f"  Lost: {lost}")
                print(f"  Gained: {gained}")
                print(f"  Stable: {len(stable)}")

                # Adjacency
                print(f"\n  T adjacency:")
                for x in range(n):
                    print(f"    {x} -> {[y for y in range(n) if A[x][y]]}")

                print(f"\n  T' adjacency:")
                for x in range(n):
                    print(f"    {x} -> {[y for y in range(n) if A2[x][y]]}")

                # Z_1 change: compute basis of Z_1 for both
                # d_1 maps arcs to vertices: d_1(u,v) = v - u
                # Z_1 = ker(d_1) in Omega_1
                arcs_T = p1[1]
                d1_mat = np.zeros((n, len(arcs_T)))
                for idx_a, (a,b) in enumerate(arcs_T):
                    d1_mat[b, idx_a] = 1
                    d1_mat[a, idx_a] = -1
                U1, S1, Vt1 = np.linalg.svd(d1_mat, full_matrices=True)
                rk1 = int(np.sum(np.abs(S1) > 1e-8))
                z1_basis_T = Vt1[rk1:].T  # columns are Z_1 basis

                arcs_T2 = p2[1]
                d1_mat2 = np.zeros((n, len(arcs_T2)))
                for idx_a, (a,b) in enumerate(arcs_T2):
                    d1_mat2[b, idx_a] = 1
                    d1_mat2[a, idx_a] = -1
                U12, S12, Vt12 = np.linalg.svd(d1_mat2, full_matrices=True)
                rk12 = int(np.sum(np.abs(S12) > 1e-8))
                z1_basis_T2 = Vt12[rk12:].T

                print(f"\n  dim(Z_1) in T: {z1_basis_T.shape[1]}")
                print(f"  dim(Z_1) in T': {z1_basis_T2.shape[1]}")

                # B_1 = im(d_2) in Z_1
                # d_2 maps TT triples to arcs
                if bdo1[2] is not None:
                    # im(d_2) in arc space
                    bd2_raw = bd1[2] @ o1[2]  # |arcs| x dim(Om2)
                    rk_bd2_T = np.linalg.matrix_rank(bd2_raw, tol=1e-8)
                else:
                    rk_bd2_T = 0

                if bdo2[2] is not None:
                    bd2_raw2 = bd2[2] @ o2[2]
                    rk_bd2_T2 = np.linalg.matrix_rank(bd2_raw2, tol=1e-8)
                else:
                    rk_bd2_T2 = 0

                print(f"  rk(d_2) in T: {rk_bd2_T} => beta_1 = {z1_basis_T.shape[1] - rk_bd2_T}")
                print(f"  rk(d_2) in T': {rk_bd2_T2} => beta_1 = {z1_basis_T2.shape[1] - rk_bd2_T2}")

                # The "new" 1-cycle: find a Z_1 element that's a boundary in T' but not T
                # Or: find a 1-cycle in T that's NOT a boundary in T, but whose image
                # IS a boundary in T'

                # Actually: beta_1(T) = dim(Z_1) - rk(d_2|T) and beta_1(T') = dim(Z_1') - rk(d_2|T')
                # dim(Z_1) = dim(Z_1') (both = (n-1)(n-2)/2)
                # beta_1 increases by 1 means rk(d_2) DECREASES by 1

                print(f"\n  KEY: rk(d_2) decreases by 1, meaning one 1-boundary is 'lost'")
                print(f"  And rk(d_3) increases by 1, meaning one 2-cycle is newly filled")

                # What 3-cycle (if any) contains the flipped arc?
                for w in range(n):
                    if w == u or w == v:
                        continue
                    # Check if {u,v,w} forms a 3-cycle in T
                    if A[u][v] and A[v][w] and A[w][u]:
                        print(f"  3-cycle in T: {u}->{v}->{w}->{u}")
                    elif A[u][w] and A[w][v] and A[v][u]:
                        print(f"  3-cycle in T: {u}->{w}->{v}->{u}")

                found = True
                break

        if found:
            break
    if found:
        break

# Now do a SECOND example: dOm2 = +1, db1 = 0 (most common non-trivial case)
print(f"\n{'='*70}")
print("EXAMPLE 2: dOm2=+1, db1=0 (dimension change, no beta_1 change)")
print("=" * 70)

found2 = False
for bits in range(total):
    A = build_adj(n, bits)
    p1, o1, bd1, bdo1, d1, r1, b1 = compute_chain_data(A, n)

    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            A2 = build_adj(n, bits2)
            p2, o2, bd2, bdo2, d2, r2, b2 = compute_chain_data(A2, n)

            dOm2 = d2[2] - d1[2]
            db1 = b2[1] - b1[1]

            if dOm2 == 1 and db1 == 0 and b1[1] == 1:
                if A[i][j] == 1:
                    u, v = i, j
                else:
                    u, v = j, i

                print(f"\nExample: bits={bits}, flip arc ({u}->{v}) to ({v}->{u})")
                print(f"  beta(T): {[b1[p] for p in range(4)]}")
                print(f"  beta(T'): {[b2[p] for p in range(4)]}")
                print(f"  dim(Om2): {d1[2]} -> {d2[2]}")
                print(f"  dim(Z2): {d1[2]-r1[2]} -> {d2[2]-r2[2]}")
                print(f"  rk(d3): {r1[3]} -> {r2[3]}")
                print(f"  rk(d2): {r1[2]} -> {r2[2]}")

                # 3-cycles
                for w in range(n):
                    if w == u or w == v:
                        continue
                    if A[u][v] and A[v][w] and A[w][u]:
                        print(f"  3-cycle killed: {u}->{v}->{w}->{u}")

                found2 = True
                break
        if found2:
            break
    if found2:
        break

# Summary of all delta combinations at n=5
print(f"\n{'='*70}")
print("COMPLETE DELTA CLASSIFICATION at n=5")
print("=" * 70)
print("(dOm2, db1, drk3, dZ2) with dZ2 = dOm2+db1 and drk3 = dZ2:")
summary = {}
for bits in range(total):
    p1, o1, bd1, bdo1, d1, r1, b1 = compute_chain_data(build_adj(n, bits), n)
    for i in range(n):
        for j in range(i+1, n):
            bits2 = flip_arc(bits, i, j, n)
            p2, o2, bd2, bdo2, d2, r2, b2 = compute_chain_data(build_adj(n, bits2), n)
            key = (d2[2]-d1[2], b2[1]-b1[1])
            if key not in summary:
                summary[key] = 0
            summary[key] += 1

for key in sorted(summary.keys()):
    dOm2, db1 = key
    dZ2 = dOm2 + db1
    print(f"  dOm2={dOm2:+d}, db1={db1:+d} => dZ2=drk3={dZ2:+d}: {summary[key]} flips")

print("\nDone.")
