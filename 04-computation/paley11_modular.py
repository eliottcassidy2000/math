#!/usr/bin/env python3
"""
paley11_modular.py — opus-2026-03-13-S71

Compute P_11 GLMY Ω dims using modular arithmetic.
Use mod small primes (97, 101) to compute ranks exactly.
"""

import numpy as np
import time
import sys

p = 11
QR = {a*a % p for a in range(1, p)}

# Adjacency
adj = [set() for _ in range(p)]
for i in range(p):
    for s in QR:
        adj[i].add((i + s) % p)

def enumerate_paths_incremental(prev_paths, adj_list, n):
    """Extend (m-1)-paths to m-paths."""
    paths = []
    for path in prev_paths:
        last = path[-1]
        visited = set(path)
        for v in adj_list[last]:
            if v not in visited:
                paths.append(path + (v,))
    return paths

def rank_mod_p(J_rows, J_cols, J_vals, nrows, ncols, mod):
    """Compute rank of sparse matrix mod a prime using Gaussian elimination.
    Much more memory efficient than dense SVD."""
    # Build dense matrix mod p
    # If too large, use sparse row reduction
    if nrows * ncols > 50_000_000:
        print(f"      Matrix too large ({nrows}x{ncols}), skipping", flush=True)
        return -1

    M = np.zeros((nrows, ncols), dtype=np.int64)
    for r, c, v in zip(J_rows, J_cols, J_vals):
        M[r, c] = (M[r, c] + v) % mod

    # Gaussian elimination mod p
    rank = 0
    for col in range(ncols):
        # Find pivot
        pivot = -1
        for row in range(rank, nrows):
            if M[row, col] % mod != 0:
                pivot = row
                break
        if pivot == -1:
            continue
        # Swap
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        # Scale
        inv = pow(int(M[rank, col]), mod - 2, mod)
        M[rank] = (M[rank] * inv) % mod
        # Eliminate
        for row in range(nrows):
            if row != rank and M[row, col] != 0:
                M[row] = (M[row] - M[row, col] * M[rank]) % mod
        rank += 1

    return rank

# Enumerate paths
print(f"P_{p}: QR = {sorted(QR)}")
print(f"\nEnumerating and computing Ω dimensions:")

MOD = 97  # Small prime for modular arithmetic

all_paths = {0: [(v,) for v in range(p)]}
omega_dims = {0: p}
print(f"  m=0: |A_0|={p}, Ω_0={p}")

for m in range(1, p):
    t0 = time.time()
    all_paths[m] = enumerate_paths_incremental(all_paths[m-1], adj, p)
    t_enum = time.time() - t0

    dim_Am = len(all_paths[m])
    print(f"\n  m={m}: |A_m|={dim_Am} (enum {t_enum:.1f}s)", flush=True)

    if m <= 1:
        omega_dims[m] = dim_Am
        print(f"         Ω_{m}={dim_Am}")
        continue

    # Build junk map
    t0 = time.time()
    am1_set = set(all_paths[m-1])

    junk = {}
    junk_count = 0
    rows_list = []
    cols_list = []
    vals_list = []

    for j, path in enumerate(all_paths[m]):
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face not in am1_set:
                if face not in junk:
                    junk[face] = junk_count
                    junk_count += 1
                rows_list.append(junk[face])
                cols_list.append(j)
                vals_list.append((-1)**i)

    t_junk = time.time() - t0
    print(f"         junk_faces={junk_count}, nnz={len(rows_list)} ({t_junk:.1f}s)", flush=True)

    if junk_count == 0:
        omega_dims[m] = dim_Am
        print(f"         Ω_{m}={dim_Am}")
        continue

    # Use the smaller dimension for rank computation
    if junk_count <= dim_Am:
        # Rank of JJ^T gives rank of J
        # JJ^T is junk_count × junk_count
        print(f"         Computing rank via JJ^T ({junk_count}x{junk_count})...", flush=True)
        t0 = time.time()

        if junk_count * junk_count > 50_000_000:
            print(f"         TOO LARGE, using sparse approach", flush=True)
            # Try computing J^TJ instead if dim_Am is smaller
            if dim_Am < junk_count and dim_Am * dim_Am <= 50_000_000:
                print(f"         Trying J^TJ ({dim_Am}x{dim_Am})...", flush=True)
                JTJ = np.zeros((dim_Am, dim_Am), dtype=np.int64)
                # Build JTJ directly: JTJ[a,b] = sum_k J[k,a]*J[k,b]
                # Group by junk face
                face_to_entries = {}
                for idx in range(len(rows_list)):
                    r = rows_list[idx]
                    c = cols_list[idx]
                    v = vals_list[idx]
                    if r not in face_to_entries:
                        face_to_entries[r] = []
                    face_to_entries[r].append((c, v))

                for r, entries in face_to_entries.items():
                    for c1, v1 in entries:
                        for c2, v2 in entries:
                            JTJ[c1, c2] = (JTJ[c1, c2] + v1 * v2) % MOD

                rank_J = np.linalg.matrix_rank(JTJ.astype(float), tol=0.5)
                t1 = time.time()
                print(f"         rank(J)={rank_J} ({t1-t0:.1f}s)", flush=True)
            else:
                print(f"         SKIPPING (both JJ^T and J^TJ too large)")
                break
        else:
            # Dense JJ^T
            JJT = np.zeros((junk_count, junk_count), dtype=np.int64)
            face_to_entries = {}
            for idx in range(len(rows_list)):
                r = rows_list[idx]
                c = cols_list[idx]
                v = vals_list[idx]
                if c not in face_to_entries:
                    face_to_entries[c] = []
                face_to_entries[c].append((r, v))

            for c, entries in face_to_entries.items():
                for r1, v1 in entries:
                    for r2, v2 in entries:
                        JJT[r1, r2] += v1 * v2

            rank_J = np.linalg.matrix_rank(JJT.astype(float), tol=0.5)
            t1 = time.time()
            print(f"         rank(J)={rank_J} ({t1-t0:.1f}s)", flush=True)
    else:
        # Rank of J^TJ gives rank of J
        print(f"         Computing rank via J^TJ ({dim_Am}x{dim_Am})...", flush=True)
        if dim_Am * dim_Am > 50_000_000:
            print(f"         TOO LARGE, skipping")
            break

        t0 = time.time()
        JTJ = np.zeros((dim_Am, dim_Am), dtype=np.int64)
        face_to_entries = {}
        for idx in range(len(rows_list)):
            r = rows_list[idx]
            c = cols_list[idx]
            v = vals_list[idx]
            if r not in face_to_entries:
                face_to_entries[r] = []
            face_to_entries[r].append((c, v))

        for r, entries in face_to_entries.items():
            for c1, v1 in entries:
                for c2, v2 in entries:
                    JTJ[c1, c2] += v1 * v2

        rank_J = np.linalg.matrix_rank(JTJ.astype(float), tol=0.5)
        t1 = time.time()
        print(f"         rank(J)={rank_J} ({t1-t0:.1f}s)", flush=True)

    omega_m = dim_Am - rank_J
    omega_dims[m] = omega_m
    div = omega_m // p if omega_m % p == 0 else f"({omega_m}/{p})"
    print(f"         Ω_{m}={omega_m} (={div}×{p})", flush=True)

    # Check memory limit
    if dim_Am > 100000:
        print(f"  Stopping: next degree will have too many paths")
        break

# Summary
print(f"\n{'='*70}")
print("SUMMARY")
print("="*70)
max_computed = max(omega_dims.keys())
om_list = [omega_dims.get(m, '?') for m in range(max_computed + 1)]
print(f"  Ω = {om_list}")
div_list = [omega_dims[m]//p if omega_dims.get(m,0)%p==0 else '?' for m in range(max_computed+1)]
print(f"  Ω/{p} = {div_list}")

print(f"\n  P_7 Ω/7 = [1, 3, 6, 9, 9, 6, 3]")
print(f"  P_11 Ω/11 (so far) = {div_list}")

# Check tail palindrome
if len(div_list) >= 3 and all(isinstance(x, int) for x in div_list):
    tail = div_list[1:]
    is_pal = tail == tail[::-1]
    print(f"  Ω_1..Ω_{max_computed}/{p} = {tail}, palindromic? {is_pal}")

print("\nDONE.")
