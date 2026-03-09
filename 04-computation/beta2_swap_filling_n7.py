#!/usr/bin/env python3
"""
beta2_swap_filling_n7.py - Test swap cycle filling at n=7 and construct explicit fillings

Key result from n=5,6: ALL swap cycles are boundaries in Omega_3.
Test this at n=7 where the Omega complex is much larger.

Also: analyze the STRUCTURE of the filling 3-chains.
What 3-paths are used? Is there a pattern?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                A[i][j] = 1
    return A


def test_swap_cycles(A, n, verbose=False):
    """Test whether all swap cycles for tournament A are boundaries.
    Returns (all_ok, num_tested, details)."""
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
    path2_idx = {p: i for i, p in enumerate(paths2)}

    all_ok = True
    num_tested = 0
    details = []

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C

        if ker_dim == 0:
            continue

        _, _, Vt = np.linalg.svd(C, full_matrices=True)
        ker_basis = Vt[rank_C:]

        for ki in range(min(ker_basis.shape[0], 5)):
            M_vec = ker_basis[ki]
            z = np.zeros(len(paths2))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12: continue
                if (a,b,v) in path2_idx and (v,a,b) in path2_idx:
                    z[path2_idx[(a,b,v)]] += coeff
                    z[path2_idx[(v,a,b)]] -= coeff

            if np.max(np.abs(z)) < 1e-12:
                continue

            num_tested += 1

            if dim_O3 > 0:
                D3 = build_full_boundary_matrix(paths3, paths2)
                D3_omega = D3 @ omega3
                w, res, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
                err = np.max(np.abs(D3_omega @ w - z))

                if err > 1e-6:
                    all_ok = False
                    details.append((v, ker_dim, err, 'FAIL'))
                elif verbose:
                    # Analyze the filling
                    w_A3 = omega3 @ w
                    nonzero_3paths = [(paths3[j], w_A3[j]) for j in range(len(paths3)) if abs(w_A3[j]) > 1e-10]
                    details.append((v, ker_dim, err, nonzero_3paths))
            else:
                if np.max(np.abs(z)) > 1e-12:
                    all_ok = False
                    details.append((v, ker_dim, 999, 'NO_O3'))

    return all_ok, num_tested, details


# ============================================================
# PART 1: n=7 random sampling
# ============================================================
print("=" * 70)
print("n=7: SWAP CYCLE FILLING TEST (500 random + Paley)")
print("=" * 70)

n = 7
t0 = time.time()
total_tested = 0
total_ok = 0
total_fail = 0

for trial in range(500):
    A = random_tournament(n)
    ok, tested, details = test_swap_cycles(A, n)
    total_tested += tested
    if ok:
        total_ok += 1
    else:
        total_fail += 1
        if total_fail <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  FAIL at trial {trial}: scores={scores}")
            for d in details:
                if d[3] == 'FAIL':
                    print(f"    v={d[0]}, ker_dim={d[1]}, err={d[2]:.2e}")

    if (trial+1) % 100 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/500 ({elapsed:.0f}s) tested={total_tested} ok={total_ok} fail={total_fail}")

# Paley T_7
A_paley = paley_tournament(7)
ok, tested, details = test_swap_cycles(A_paley, n, verbose=True)
total_tested += tested
elapsed = time.time() - t0

print(f"\nn=7 ({elapsed:.0f}s): {total_tested} swap cycles tested")
print(f"  All OK: {total_ok}/500 random + Paley={'OK' if ok else 'FAIL'}")
print(f"  Failures: {total_fail}")

if ok and details:
    print(f"\n  Paley T_7 swap cycle details:")
    for d in details[:3]:
        v, ker_dim, err, filling = d
        if isinstance(filling, list):
            print(f"    v={v}: ker_dim={ker_dim}, err={err:.2e}, filling_size={len(filling)}")
            for path, coeff in filling[:5]:
                print(f"      {path}: {coeff:+.4f}")


# ============================================================
# PART 2: Analyze structure of fillings at n=5
# ============================================================
print(f"\n{'='*70}")
print("n=5: FILLING STRUCTURE ANALYSIS")
print("=" * 70)

n = 5


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


filling_types = defaultdict(int)  # {(v_positions): count}

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)

    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C

        if ker_dim == 0:
            continue

        _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
        ker_basis = Vt[rank_C:]

        paths2 = enumerate_allowed_paths(A, n, 2)
        paths3 = enumerate_allowed_paths(A, n, 3)
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
        path2_idx = {p: i for i, p in enumerate(paths2)}

        if dim_O3 == 0:
            continue

        D3 = build_full_boundary_matrix(paths3, paths2)
        D3_omega = D3 @ omega3

        for ki in range(ker_basis.shape[0]):
            M_vec = ker_basis[ki]
            z = np.zeros(len(paths2))
            for j, (a, b) in enumerate(arcs_PQ):
                coeff = M_vec[j]
                if abs(coeff) < 1e-12: continue
                if (a,b,v) in path2_idx and (v,a,b) in path2_idx:
                    z[path2_idx[(a,b,v)]] += coeff
                    z[path2_idx[(v,a,b)]] -= coeff

            if np.max(np.abs(z)) < 1e-12:
                continue

            w, res, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
            err = np.max(np.abs(D3_omega @ w - z))

            if err < 1e-6:
                w_A3 = omega3 @ w
                # Classify 3-paths by position of v
                v_positions = []
                for j, path in enumerate(paths3):
                    if abs(w_A3[j]) > 1e-10:
                        if path[0] == v: v_positions.append(0)
                        elif path[1] == v: v_positions.append(1)
                        elif path[2] == v: v_positions.append(2)
                        elif path[3] == v: v_positions.append(3)
                        else: v_positions.append(-1)  # v not in path

                filling_types[tuple(sorted(set(v_positions)))] += 1

print(f"\nFilling structure (position of v in 3-paths):")
for positions, count in sorted(filling_types.items(), key=lambda x: -x[1]):
    pos_labels = {0: 'front', 1: 'pos1', 2: 'pos2', 3: 'back', -1: 'absent'}
    labels = [pos_labels.get(p, str(p)) for p in positions]
    print(f"  v at positions {positions} ({', '.join(labels)}): {count}")


# ============================================================
# PART 3: Explicit filling for one example
# ============================================================
print(f"\n{'='*70}")
print("EXPLICIT FILLING EXAMPLE")
print("=" * 70)

# Find a specific swap cycle at n=5 and show the full filling
n = 5
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    for v in range(n):
        P = [a for a in range(n) if a != v and A[v][a] == 1]
        Q = [b for b in range(n) if b != v and A[b][v] == 1]
        arcs_PQ = [(a, b) for a in P for b in Q if A[a][b] == 1]

        if len(arcs_PQ) < 2:
            continue

        m = len(arcs_PQ)
        rows = []
        for a in P:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if a2 == a: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)
        for b in Q:
            row = [0] * m
            for j, (a2, b2) in enumerate(arcs_PQ):
                if b2 == b: row[j] = 1
            if any(r != 0 for r in row):
                rows.append(row)

        if not rows:
            continue

        C_mat = np.array(rows, dtype=float)
        Sc = np.linalg.svd(C_mat, compute_uv=False)
        rank_C = sum(s > 1e-8 for s in Sc)
        ker_dim = m - rank_C

        if ker_dim == 0:
            continue

        _, _, Vt = np.linalg.svd(C_mat, full_matrices=True)
        M_vec = Vt[rank_C]

        paths2 = enumerate_allowed_paths(A, n, 2)
        paths3 = enumerate_allowed_paths(A, n, 3)
        omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
        dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0
        path2_idx = {p: i for i, p in enumerate(paths2)}

        z = np.zeros(len(paths2))
        for j, (a, b) in enumerate(arcs_PQ):
            coeff = M_vec[j]
            if abs(coeff) < 1e-12: continue
            if (a,b,v) in path2_idx and (v,a,b) in path2_idx:
                z[path2_idx[(a,b,v)]] += coeff
                z[path2_idx[(v,a,b)]] -= coeff

        if np.max(np.abs(z)) < 1e-12:
            continue

        if dim_O3 == 0:
            continue

        D3 = build_full_boundary_matrix(paths3, paths2)
        D3_omega = D3 @ omega3
        w, res, _, _ = np.linalg.lstsq(D3_omega, z, rcond=None)
        err = np.max(np.abs(D3_omega @ w - z))

        if err < 1e-6:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"T#{bits} scores={scores}, v={v}")
            print(f"  P={P}, Q={Q}")
            print(f"  Arcs P->Q: {arcs_PQ}")
            print(f"  Swap cycle M coefficients:")
            for j, (a, b) in enumerate(arcs_PQ):
                if abs(M_vec[j]) > 1e-10:
                    print(f"    M[{a},{b}] = {M_vec[j]:+.6f}")

            print(f"\n  Swap cycle z (A_2 coords):")
            for j, path in enumerate(paths2):
                if abs(z[j]) > 1e-10:
                    print(f"    z[{path}] = {z[j]:+.6f}")

            print(f"\n  Filling w (A_3 coords), err={err:.2e}:")
            w_A3 = omega3 @ w
            for j, path in enumerate(paths3):
                if abs(w_A3[j]) > 1e-10:
                    v_pos = path.index(v) if v in path else -1
                    print(f"    w[{path}] = {w_A3[j]:+.6f} (v at pos {v_pos})")

            # Verify: d_3(w) = z
            check = D3 @ w_A3
            print(f"\n  Verification d_3(w) = z: max_err = {np.max(np.abs(check - z)):.2e}")

            # Check: is the SWAP cycle actually in Omega_2?
            # It should be by construction, but let's verify
            paths1 = enumerate_allowed_paths(A, n, 1)
            omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
            # Project z onto Omega_2: z should be in the column space of omega2
            proj_z, _, _, _ = np.linalg.lstsq(omega2, z, rcond=None)
            z_projected = omega2 @ proj_z
            omega_err = np.max(np.abs(z - z_projected))
            print(f"  z in Omega_2? err = {omega_err:.2e}")

            break
    else:
        continue
    break


print("\n\nDone.")
