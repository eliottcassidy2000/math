#!/usr/bin/env python3
"""
beta2_full_chain.py - Full Omega chain complex analysis

Compute dim(Omega_p) for ALL p and ALL tournaments at small n.
Look for patterns, Euler characteristic, and structural constraints
that explain why H_2 = 0 universally.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from itertools import permutations, combinations
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def full_omega_dims(A, n, max_p=None):
    """Compute dim(Omega_p) for all p up to max_p."""
    if max_p is None:
        max_p = n - 1

    allowed = {}
    for p in range(-1, max_p + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    omega_dims = []
    for p in range(max_p + 2):
        omega = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        dim = omega.shape[1] if omega.ndim == 2 and omega.shape[1] > 0 else 0
        omega_dims.append(dim)

    return omega_dims


def full_betti(A, n, max_p=None):
    """Compute ALL Betti numbers up to max_p."""
    if max_p is None:
        max_p = n - 1
    return path_betti_numbers(A, n, max_dim=max_p)


def count_3cycles(A, n):
    c3 = 0
    for triple in combinations(range(n), 3):
        i,j,k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
            c3 += 1
    return c3


# ============================================================
# PART 1: Full chain complex at n=4,5
# ============================================================
print("=" * 70)
print("FULL OMEGA CHAIN COMPLEX DIMENSIONS")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    score_omega = defaultdict(list)
    score_betti = defaultdict(list)

    for bits in range(1 << (n*(n-1)//2)):
        A = build_adj(n, bits)
        scores = tuple(sorted([sum(row) for row in A]))
        omega_dims = full_omega_dims(A, n)
        betti = full_betti(A, n)
        score_omega[scores].append(omega_dims)
        score_betti[scores].append(betti)

    for scores in sorted(score_omega.keys()):
        items = score_omega[scores]
        bettis = score_betti[scores]
        cnt = len(items)

        # Check if all have same omega dims
        all_same_omega = all(x == items[0] for x in items)
        all_same_betti = all(x == bettis[0] for x in bettis)

        print(f"\n  Score {scores} ({cnt} tournaments):")
        if all_same_omega:
            print(f"    Omega dims: {items[0]}")
        else:
            # Show distribution
            omega_dist = Counter(tuple(x) for x in items)
            for od, c in sorted(omega_dist.items()):
                print(f"    Omega dims: {list(od)} x{c}")
        if all_same_betti:
            print(f"    Betti: {bettis[0]}")
        else:
            betti_dist = Counter(tuple(x) for x in bettis)
            for b, c in sorted(betti_dist.items()):
                print(f"    Betti: {list(b)} x{c}")

        # Euler char
        for item in items[:1]:
            chi_omega = sum((-1)**p * item[p] for p in range(len(item)))
            chi_betti = sum((-1)**p * bettis[0][p] for p in range(len(bettis[0])))
            print(f"    Euler(omega): {chi_omega}, Euler(betti): {chi_betti}")


# ============================================================
# PART 2: Explicit 2-cycles and their fillings at n=5
# ============================================================
print(f"\n{'='*70}")
print("EXPLICIT 2-CYCLES AND 3-CHAIN FILLINGS (n=5)")
print("=" * 70)

n = 5

# Take one regular tournament
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = sorted([sum(row) for row in A])
    if scores == [2,2,2,2,2]:
        break

print(f"\nRegular tournament T#{bits}")
print(f"Adjacency:")
for i in range(n):
    row = ' '.join(str(A[i][j]) for j in range(n))
    print(f"  {i}: {row}")

# Compute Omega_2, Omega_3 bases
paths2 = enumerate_allowed_paths(A, n, 2)
paths3 = enumerate_allowed_paths(A, n, 3)
paths1 = enumerate_allowed_paths(A, n, 1)

omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
print(f"\ndim(Omega_2) = {omega2.shape[1]}")
print(f"dim(Omega_3) = {omega3.shape[1]}")

# Compute boundary d2: Omega_2 -> A_1
path1_idx = {p: i for i, p in enumerate(paths1)}
D2_full = np.zeros((len(paths1), len(paths2)))
for j, (a,b,c) in enumerate(paths2):
    if (b,c) in path1_idx: D2_full[path1_idx[(b,c)], j] += 1
    if (a,c) in path1_idx: D2_full[path1_idx[(a,c)], j] -= 1
    if (a,b) in path1_idx: D2_full[path1_idx[(a,b)], j] += 1

D2_omega = D2_full @ omega2
U2, S2, V2t = np.linalg.svd(D2_omega, full_matrices=True)
rank_d2 = sum(s > 1e-8 for s in S2)
ker_d2_dim = omega2.shape[1] - rank_d2
print(f"rank(d_2) = {rank_d2}, ker(d_2) = {ker_d2_dim}")

# Get kernel basis of d_2
# Kernel in Omega_2 coordinates
ker_d2_basis = V2t[rank_d2:].T  # columns are kernel vectors in Omega_2 coordinates

# Convert to A_2 coordinates
ker_d2_A2 = omega2 @ ker_d2_basis  # columns in A_2 coordinates

print(f"\n2-cycle generators (in A_2 coordinates):")
for k in range(min(ker_d2_dim, 3)):
    v = ker_d2_A2[:, k]
    nonzero = [(paths2[i], round(v[i], 4)) for i in range(len(v)) if abs(v[i]) > 1e-6]
    print(f"  z_{k}: {len(nonzero)} nonzero coefficients")
    for path, coeff in nonzero[:10]:
        print(f"    {path}: {coeff}")
    if len(nonzero) > 10:
        print(f"    ... and {len(nonzero)-10} more")

# Boundary d3: Omega_3 -> A_2
D3_full = build_full_boundary_matrix(paths3, paths2)
D3_omega = D3_full @ omega3
im_d3_dim = sum(s > 1e-8 for s in np.linalg.svd(D3_omega, compute_uv=False))
print(f"\nim(d_3) = {im_d3_dim}")
print(f"beta_2 = ker - im = {ker_d2_dim} - {im_d3_dim} = {ker_d2_dim - im_d3_dim}")


# ============================================================
# PART 3: CONE CONSTRUCTION TEST
# ============================================================
print(f"\n{'='*70}")
print("CONE CONSTRUCTION: CAN WE FILL 2-CYCLES?")
print("=" * 70)

# For a 2-cycle z, try to construct a 3-chain w with d_3(w) = z
# using the cone construction with vertex v.

# For each 2-path (a,b,c) in z, and fixed cone vertex v:
# Try all possible 3-paths that "cone" (a,b,c) with v:
# (v,a,b,c), (a,v,b,c), (a,b,v,c), (a,b,c,v)
# but only those that are allowed

path3_idx = {p: i for i, p in enumerate(paths3)}

for cone_v in range(n):
    print(f"\nCone vertex v={cone_v}:")

    # For each 2-path, find which 3-paths involve cone_v
    cone_map = {}  # 2-path -> list of (3-path, sign, position)
    for i, p3 in enumerate(paths3):
        if cone_v not in p3:
            continue
        pos = list(p3).index(cone_v)
        # The face at position pos (deleting cone_v) gives back a 2-path
        face = p3[:pos] + p3[pos+1:]
        sign = (-1)**pos
        if face in {p: None for p in paths2}:  # check if face is a 2-path
            face_idx = next(j for j, p in enumerate(paths2) if p == face)
            if face not in cone_map:
                cone_map[face] = []
            cone_map[face].append((p3, sign, pos))

    covered = len(cone_map)
    total = len(paths2)
    print(f"  2-paths coverable by cone: {covered}/{total}")

    # Check: for the first 2-cycle, can we write it as d_3(cone)?
    if ker_d2_dim > 0:
        z = ker_d2_A2[:, 0]  # first 2-cycle

        # Check which 2-paths in z have nonzero coefficient
        nonzero_paths = [(i, paths2[i], z[i]) for i in range(len(z)) if abs(z[i]) > 1e-6]
        in_cone = sum(1 for i, p, c in nonzero_paths if p in cone_map)
        print(f"  2-cycle z_0: {len(nonzero_paths)} nonzero, {in_cone} in cone domain")


# ============================================================
# PART 4: Is dim(Omega_p) determined by score sequence?
# ============================================================
print(f"\n{'='*70}")
print("IS dim(Omega_p) DETERMINED BY SCORE SEQUENCE?")
print("=" * 70)

n = 5
score_to_omega = defaultdict(set)
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    omega_dims = tuple(full_omega_dims(A, n))
    score_to_omega[scores].add(omega_dims)

all_determined = True
for scores in sorted(score_to_omega.keys()):
    vals = score_to_omega[scores]
    determined = len(vals) == 1
    if not determined:
        all_determined = False
    print(f"  Score {scores}: {len(vals)} distinct Omega profiles {'(UNIQUE)' if determined else '(VARIES!)'}")
    if not determined:
        for v in sorted(vals):
            print(f"    {list(v)}")

print(f"\n  Omega profile determined by score? {all_determined}")


# ============================================================
# PART 5: Formula for Euler characteristic
# ============================================================
print(f"\n{'='*70}")
print("EULER CHARACTERISTIC chi = sum (-1)^p dim(Omega_p)")
print("=" * 70)

for n in [3, 4, 5]:
    chi_vals = Counter()
    for bits in range(1 << (n*(n-1)//2)):
        A = build_adj(n, bits)
        omega_dims = full_omega_dims(A, n)
        chi = sum((-1)**p * omega_dims[p] for p in range(len(omega_dims)))
        chi_vals[chi] += 1

    print(f"\nn={n}: chi distribution: {dict(sorted(chi_vals.items()))}")
    print(f"  chi = 1 always? {'YES' if len(chi_vals)==1 and 1 in chi_vals else 'NO'}")


# ============================================================
# PART 6: beta_1 relationship to Omega_2 kernel
# ============================================================
print(f"\n{'='*70}")
print("RELATIONSHIP: beta_1, c3, ker(d_2), im(d_3)")
print("=" * 70)

n = 5
data = []
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = count_3cycles(A, n)
    betti = path_betti_numbers(A, n, max_dim=3)

    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_omega2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_omega3 = omega3.shape[1] if omega3.ndim == 2 else 0

    # Boundary ranks
    if dim_omega2 > 0:
        D2 = np.zeros((len(paths1), len(paths2)))
        p1_idx = {p: i for i, p in enumerate(paths1)}
        for j, (a,b,c) in enumerate(paths2):
            if (b,c) in p1_idx: D2[p1_idx[(b,c)], j] += 1
            if (a,c) in p1_idx: D2[p1_idx[(a,c)], j] -= 1
            if (a,b) in p1_idx: D2[p1_idx[(a,b)], j] += 1
        D2_o = D2 @ omega2
        rank_d2 = sum(s > 1e-8 for s in np.linalg.svd(D2_o, compute_uv=False))
    else:
        rank_d2 = 0

    ker_d2 = dim_omega2 - rank_d2

    data.append({
        'c3': c3, 'b1': betti[1], 'b2': betti[2] if len(betti)>2 else 0,
        'dim_omega2': dim_omega2, 'dim_omega3': dim_omega3,
        'rank_d2': rank_d2, 'ker_d2': ker_d2,
    })

# Group by (c3, beta_1)
from collections import defaultdict
grouped = defaultdict(list)
for d in data:
    grouped[(d['c3'], d['b1'])].append(d)

print(f"\nn=5: (c3, beta_1) -> (dim_Omega2, ker_d2, dim_Omega3)")
for key in sorted(grouped.keys()):
    items = grouped[key]
    omega2_vals = set(d['dim_omega2'] for d in items)
    ker_vals = set(d['ker_d2'] for d in items)
    omega3_vals = set(d['dim_omega3'] for d in items)
    print(f"  c3={key[0]}, b1={key[1]}: ({sorted(omega2_vals)}, {sorted(ker_vals)}, {sorted(omega3_vals)}) [{len(items)} tours]")


# ============================================================
# PART 7: KEY TEST - Is ker(d_2) = beta_1 + C(n,4)?
# ============================================================
print(f"\n{'='*70}")
print("FORMULA TEST: ker(d_2) = f(beta_1, n, c3)?")
print("=" * 70)

for d in data[:5]:
    print(f"  c3={d['c3']}, b1={d['b1']}: ker_d2={d['ker_d2']}, C(5,4)={5-1}, b1+C(n-1,3)={d['b1']+4}")

# Check: is ker(d_2) related to c3 or beta_1?
# For transitive: ker=4 = C(5,4)
# For regular: ker=5 = C(5,4) + 1 = C(5,4) + beta_1

formula_match = True
for d in data:
    expected = 4 + d['b1']  # = C(5,4) + beta_1
    if d['ker_d2'] != expected:
        formula_match = False
        break

print(f"\n  ker(d_2) = C(n,4) + beta_1 = 4 + beta_1? {formula_match}")

if not formula_match:
    print("  Checking other formulas...")
    # Try ker_d2 = c3 + something
    for d in data[:20]:
        print(f"    c3={d['c3']}, b1={d['b1']}, ker={d['ker_d2']}")


print("\n\nDone.")
