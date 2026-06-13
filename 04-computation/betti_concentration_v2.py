#!/usr/bin/env python3
"""
BETTI CONCENTRATION v2 — Targeted tests
- n=6: only need max_dim=3 (β₃ first appears at n=7)
- n=7: max_dim=5 sampled (β₃ and β₄)
- n=8: max_dim=5 sampled (β₃, β₄, β₅)

Conjecture: At most one β_p (p≥1) is nonzero.

opus-2026-03-15-S72d
"""
import numpy as np
from collections import Counter
import random, time

# Inline path homology functions (no validation)
def enumerate_allowed_paths(A, n, p):
    if p < 0: return []
    if p == 0: return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1: adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(len(path))]

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx = {p: i for i, p in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx: M[idx[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    d = len(allowed_p)
    if d == 0: return np.zeros((0, 0))
    if p == 0: return np.eye(d)
    allowed_set = set(allowed_pm1)
    na_faces = {}; na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_set:
                if face not in na_faces:
                    na_faces[face] = na_count; na_count += 1
    if na_count == 0: return np.eye(d)
    P = np.zeros((na_count, d))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in na_faces: P[na_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    ns = Vt[rank:].T
    return ns if ns.shape[1] > 0 else np.zeros((d, 0))

def path_betti(A, n, max_dim):
    allowed = {p: ([] if p < 0 else enumerate_allowed_paths(A, n, p))
               for p in range(-1, max_dim + 2)}
    omega = {p: compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
             for p in range(max_dim + 2)}
    betti = []
    for p in range(max_dim + 1):
        dop = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dop == 0: betti.append(0); continue
        bd = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bdo = bd @ omega[p]
        rp = sum(s > 1e-8 for s in np.linalg.svd(bdo, compute_uv=False)) if min(bdo.shape) > 0 else 0
        kd = dop - rp
        dop1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dop1 > 0:
            bd1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bdo1 = bd1 @ omega[p+1]
            im = sum(s > 1e-8 for s in np.linalg.svd(bdo1, compute_uv=False)) if min(bdo1.shape) > 0 else 0
        else: im = 0
        betti.append(max(0, kd - im))
    return betti

def tourn(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

print("=" * 70)
print("BETTI CONCENTRATION v2")
print("=" * 70)

# n=6 EXHAUSTIVE with max_dim=3
print("\n--- n=6 EXHAUSTIVE (max_dim=3) ---")
n = 6; E = n*(n-1)//2
bd = Counter(); multi = []; t0 = time.time()
for bits in range(1 << E):
    A = tourn(n, bits)
    b = path_betti(A, n, 3)
    bt = tuple(b)
    bd[bt] += 1
    nz = [p for p in range(1, len(b)) if b[p] > 0]
    if len(nz) > 1: multi.append((bits, b))
    if (bits+1) % 10000 == 0:
        print(f"  {bits+1}/{1<<E} ({100*(bits+1)/(1<<E):.0f}%), t={time.time()-t0:.0f}s", flush=True)
print(f"  Done in {time.time()-t0:.1f}s")
for bt in sorted(bd.keys()): print(f"    β={list(bt)}: {bd[bt]}")
print(f"  Multi-nonzero: {len(multi)}")

# n=7 SAMPLED with max_dim=5
print("\n--- n=7 SAMPLED (3000, max_dim=5) ---")
n = 7; E = n*(n-1)//2
random.seed(42)
bd = Counter(); multi = []; nzd = Counter(); t0 = time.time()
for trial in range(3000):
    bits = random.randint(0, (1 << E) - 1)
    A = tourn(n, bits)
    b = path_betti(A, n, 5)
    bt = tuple(b)
    bd[bt] += 1
    nz = tuple(p for p in range(1, len(b)) if b[p] > 0)
    nzd[nz] += 1
    if len(nz) > 1:
        multi.append((bits, b))
    if (trial+1) % 500 == 0:
        print(f"  {trial+1}/3000, t={time.time()-t0:.0f}s, multi={len(multi)}", flush=True)
print(f"  Done in {time.time()-t0:.1f}s")
print(f"  Betti types:")
for bt in sorted(bd.keys()): print(f"    β={list(bt)}: {bd[bt]}")
print(f"  Nonzero dims: {dict(sorted(nzd.items()))}")
print(f"  Multi-nonzero: {len(multi)}")
if multi:
    for bits, b in multi[:10]:
        print(f"    VIOLATION: β={b}")

# n=8 SAMPLED with max_dim=5
print("\n--- n=8 SAMPLED (500, max_dim=5) ---")
n = 8; E = n*(n-1)//2
random.seed(99)
bd = Counter(); multi = []; nzd = Counter(); t0 = time.time()
for trial in range(500):
    bits = random.randint(0, (1 << E) - 1)
    A = tourn(n, bits)
    try:
        b = path_betti(A, n, 5)
        bt = tuple(b)
        bd[bt] += 1
        nz = tuple(p for p in range(1, len(b)) if b[p] > 0)
        nzd[nz] += 1
        if len(nz) > 1:
            multi.append((bits, b))
    except Exception as e:
        print(f"  Error: {e}")
    if (trial+1) % 100 == 0:
        print(f"  {trial+1}/500, t={time.time()-t0:.0f}s, multi={len(multi)}", flush=True)
print(f"  Done in {time.time()-t0:.1f}s")
print(f"  Betti types:")
for bt in sorted(bd.keys()): print(f"    β={list(bt)}: {bd[bt]}")
print(f"  Nonzero dims: {dict(sorted(nzd.items()))}")
print(f"  Multi-nonzero: {len(multi)}")
if multi:
    for bits, b in multi[:10]:
        print(f"    VIOLATION: β={b}")

# Euler characteristic for all collected
print("\n" + "=" * 70)
print("EULER CHARACTERISTIC ANALYSIS")
print("=" * 70)
for bt, cnt in sorted(bd.items()):
    chi = sum((-1)**p * v for p, v in enumerate(bt))
    print(f"  β={list(bt)}: χ={chi}, count={cnt}")

print("\nDone.")
