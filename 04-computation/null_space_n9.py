#!/usr/bin/env python3
"""
Null space of the forward-edge distribution map at n=9.

Verifies THM-065 predictions:
  - 7 OCF invariants: t3(f=6), t5(f=4), t7(f=2), t9(f=0),
                       bc33(f=4), bc35(f=2), a3(f=2)
  - 4 distinct f-values: {6, 4, 2, 0}
  - Null space dimension: 7 - 4 = 3
  - f=4 null vector: (-2)*t5 + (1)*bc33
  - f=2 null vectors: (-2)*t7 + (1)*bc35 and (-4)*t7 + (1)*a3

CRITICAL: The OCF invariants count DIRECTED CYCLES as nodes of Omega(T),
not vertex sets. A 5-vertex set can have up to 3 directed 5-cycles,
each becoming a separate node in Omega. The independence sets are in
this directed-cycle conflict graph.

kind-pasteur-2026-03-07
"""

import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from itertools import combinations, permutations
from collections import defaultdict
from fractions import Fraction
import random
import numpy as np

# =====================================================================
# Tournament generation
# =====================================================================

def random_tournament(n, seed):
    """Generate random n-tournament."""
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

# =====================================================================
# Forward-edge distribution
# =====================================================================

def count_forward_edge_dist(A, n):
    """Count a_k = number of Hamiltonian paths with exactly k forward edges."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            for k in range(n):
                c = dp.get((mask, v, k), 0)
                if c == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    new_k = k + A[v][u]
                    key = (mask | (1 << u), u, new_k)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = [0] * n
    for v in range(n):
        for k in range(n):
            dist[k] += dp.get((full, v, k), 0)
    return dist

# =====================================================================
# Directed cycle enumeration (nodes of Omega)
# =====================================================================

def find_all_directed_cycles(A, n):
    """Find ALL directed odd cycles in T. Each cycle is represented as a
    tuple of vertices in cycle order (v0 -> v1 -> ... -> v_{k-1} -> v0).

    We canonicalize: smallest vertex first, second vertex < last vertex.
    Returns list of (cycle_tuple, vertex_frozenset) pairs.
    """
    cycles = []

    for cycle_len in [3, 5, 7, 9]:
        if cycle_len > n:
            break
        for verts in combinations(range(n), cycle_len):
            # Find all directed Hamiltonian cycles on this vertex set
            # DP: dp[mask][v] = number of paths from verts[0] to v using mask
            sub = [[A[verts[i]][verts[j]] for j in range(cycle_len)] for i in range(cycle_len)]

            # Enumerate actual cycles, not just count them
            # Use DP to find all Hamiltonian paths starting at vertex 0
            # dp[mask][v] = list of paths (as tuples)
            # Too expensive to store all paths for large cycles.
            # Instead: count directed Hamiltonian cycles by DP, and for
            # the conflict graph we only need the vertex set (two cycles
            # on the same vertex set always conflict).
            #
            # WAIT: Two directed cycles on the SAME vertex set always share
            # vertices, so they always conflict in Omega. For the independence
            # polynomial, independent sets can only contain cycles on DISJOINT
            # vertex sets. So alpha_k counts collections of k pairwise
            # vertex-disjoint directed odd cycles.
            #
            # But t3 should count the number of DIRECTED 3-cycles, not just
            # vertex sets. Each cyclic triple has exactly 1 directed 3-cycle.
            # A 5-vertex set can have 0, 1, 2, or 3 directed 5-cycles.
            # For the independence polynomial I(Omega, x), alpha_1 counts
            # INDIVIDUAL directed cycles (nodes of Omega), not vertex sets.

            dp_count = [[0] * cycle_len for _ in range(1 << cycle_len)]
            dp_count[1][0] = 1
            for m in range(1, 1 << cycle_len):
                for v in range(cycle_len):
                    if not (m & (1 << v)) or dp_count[m][v] == 0:
                        continue
                    for u in range(cycle_len):
                        if m & (1 << u):
                            continue
                        if sub[v][u]:
                            dp_count[m | (1 << u)][u] += dp_count[m][v]
            full = (1 << cycle_len) - 1
            num_dir_cycles = sum(dp_count[full][v] for v in range(1, cycle_len) if sub[v][0])

            if num_dir_cycles > 0:
                vset = frozenset(verts)
                cycles.append((cycle_len, vset, num_dir_cycles))

    return cycles

def compute_ocf_invariants_correct(A, n):
    """Compute OCF invariants correctly using the conflict graph Omega(T).

    Omega(T) has one node per directed odd cycle. Two nodes are adjacent
    iff their cycles share a vertex. Since cycles on the same vertex set
    always share vertices, the independence polynomial depends only on:

    For alpha_1: sum of (number of directed k-cycles) over all vertex sets
    For alpha_2: sum over pairs of disjoint vertex sets of (product of cycle counts)
    For alpha_3: sum over triples of pairwise disjoint vertex sets of (product of cycle counts)

    This is because the conflict graph restricted to cycles on the same vertex
    set is a clique, so any independent set picks at most one cycle per vertex set.
    And for disjoint vertex sets, any pair of cycles (one from each) is independent.

    So: alpha_k = sum over k-tuples of pairwise disjoint cycle-bearing vertex sets
                  of product of (number of directed cycles on each vertex set).
    """
    all_cycles = find_all_directed_cycles(A, n)

    # Group by cycle length and vertex set
    # For each vertex set, store (cycle_length, num_directed_cycles)
    by_vset = {}  # vset -> (cycle_len, count)
    for clen, vset, count in all_cycles:
        by_vset[vset] = (clen, count)

    # Separate by cycle length
    c3_data = {}  # vset -> num_directed_cycles for 3-cycles
    c5_data = {}
    c7_data = {}
    c9_data = {}

    for vset, (clen, count) in by_vset.items():
        if clen == 3: c3_data[vset] = count
        elif clen == 5: c5_data[vset] = count
        elif clen == 7: c7_data[vset] = count
        elif clen == 9: c9_data[vset] = count

    # t_k = total number of directed k-cycles = sum of counts
    t3 = sum(c3_data.values())
    t5 = sum(c5_data.values())
    t7 = sum(c7_data.values())
    t9 = sum(c9_data.values())

    # For alpha_2 contributions:
    # bc33 = sum over disjoint 3-vset pairs of (count1 * count2)
    c3_items = list(c3_data.items())
    bc33 = 0
    for i in range(len(c3_items)):
        for j in range(i+1, len(c3_items)):
            if c3_items[i][0].isdisjoint(c3_items[j][0]):
                bc33 += c3_items[i][1] * c3_items[j][1]

    # bc35 = sum over disjoint (3-vset, 5-vset) pairs of (count3 * count5)
    c5_items = list(c5_data.items())
    bc35 = 0
    for vs3, cnt3 in c3_items:
        for vs5, cnt5 in c5_items:
            if vs3.isdisjoint(vs5):
                bc35 += cnt3 * cnt5

    # a3 = sum over triples of pairwise disjoint 3-vsets of (count1 * count2 * count3)
    a3 = 0
    for i in range(len(c3_items)):
        for j in range(i+1, len(c3_items)):
            if not c3_items[i][0].isdisjoint(c3_items[j][0]):
                continue
            for k in range(j+1, len(c3_items)):
                if c3_items[i][0].isdisjoint(c3_items[k][0]) and c3_items[j][0].isdisjoint(c3_items[k][0]):
                    a3 += c3_items[i][1] * c3_items[j][1] * c3_items[k][1]

    return t3, t5, t7, t9, bc33, bc35, a3

# =====================================================================
# MAIN COMPUTATION
# =====================================================================

n = 9
NUM_TOURNAMENTS = 250
inv_names = ['t3', 't5', 't7', 't9', 'bc33', 'bc35', 'a3']
f_values  = [  6,    4,    2,    0,     4,      2,     2  ]
p_values  = [  1,    1,    1,    1,     2,      2,     3  ]

print("=" * 72)
print(f"NULL SPACE OF FORWARD-EDGE DISTRIBUTION MAP AT n={n}")
print(f"Verifying THM-065 (f-grouping) with {NUM_TOURNAMENTS} random tournaments")
print("=" * 72)

print(f"\nOCF invariants and their f-values:")
for name, f, p in zip(inv_names, f_values, p_values):
    print(f"  {name:>5s}: f={f}, parts={p}")

print(f"\nf-groups:")
print(f"  f=6: t3 alone")
print(f"  f=4: t5, bc33        => 1 null vector")
print(f"  f=2: t7, bc35, a3    => 2 null vectors")
print(f"  f=0: t9 alone")
print(f"  Predicted null dim = 1 + 2 = 3")

# =====================================================================
# STEP 1: Generate tournaments and compute invariants + distributions
# =====================================================================

print(f"\n{'=' * 72}")
print("STEP 1: Computing invariants and forward-edge distributions")
print("=" * 72)

X = []  # [1, t3, t5, t7, t9, bc33, bc35, a3]
Y = []  # [a_0, ..., a_8]

for trial in range(NUM_TOURNAMENTS):
    if trial % 50 == 0:
        print(f"  Processing tournament {trial}/{NUM_TOURNAMENTS}...")
    A = random_tournament(n, seed=99000 + trial)
    t3, t5, t7, t9, bc33, bc35, a3 = compute_ocf_invariants_correct(A, n)
    dist = count_forward_edge_dist(A, n)

    X.append([1, t3, t5, t7, t9, bc33, bc35, a3])
    Y.append(dist)

X = np.array(X, dtype=float)
Y = np.array(Y, dtype=float)
print(f"  Done. X shape: {X.shape}, Y shape: {Y.shape}")

# =====================================================================
# STEP 2: Verify OCF: H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3
# where alpha_1 = t3+t5+t7+t9, alpha_2 = bc33+bc35, alpha_3 = a3
# =====================================================================

print(f"\n{'=' * 72}")
print("STEP 2: OCF verification (H(T) = I(Omega(T), 2))")
print("=" * 72)

ocf_errors = 0
for trial in range(NUM_TOURNAMENTS):
    inv = X[trial, 1:]
    H_ocf = 1 + 2*(inv[0]+inv[1]+inv[2]+inv[3]) + 4*(inv[4]+inv[5]) + 8*inv[6]
    H_actual = Y[trial, n-1]
    if abs(H_ocf - H_actual) > 0.5:
        ocf_errors += 1
        if ocf_errors <= 3:
            print(f"  MISMATCH trial {trial}: OCF={H_ocf:.0f}, actual={H_actual:.0f}, "
                  f"t3={int(inv[0])}, t5={int(inv[1])}, t7={int(inv[2])}, t9={int(inv[3])}, "
                  f"bc33={int(inv[4])}, bc35={int(inv[5])}, a3={int(inv[6])}")

if ocf_errors == 0:
    print(f"  OCF verified for all {NUM_TOURNAMENTS} tournaments: PASS")
    for trial in range(3):
        inv = X[trial, 1:]
        H = int(Y[trial, n-1])
        print(f"    T{trial}: t3={int(inv[0])}, t5={int(inv[1])}, t7={int(inv[2])}, "
              f"t9={int(inv[3])}, bc33={int(inv[4])}, bc35={int(inv[5])}, "
              f"a3={int(inv[6])}, H={H}")
else:
    print(f"  OCF FAILED for {ocf_errors}/{NUM_TOURNAMENTS} tournaments!")

# =====================================================================
# STEP 3: Regression a_k = c0 + sum c_I * I(T)
# =====================================================================

print(f"\n{'=' * 72}")
print("STEP 3: Linear regression a_k = c0 + sum c_I * I(T)")
print("=" * 72)

all_names = ['const'] + inv_names

header = f"{'k':>3s}"
for name in all_names:
    header += f" {name:>10s}"
header += f" {'max_err':>10s}"
print(header)

C_full = []
for k in range(n):
    coeffs, _, _, _ = np.linalg.lstsq(X, Y[:, k], rcond=None)
    pred = X @ coeffs
    max_err = np.max(np.abs(Y[:, k] - pred))
    C_full.append(coeffs)

    row = f"{k:3d}"
    for val in coeffs:
        if abs(val) < 0.001:
            row += f" {'0':>10s}"
        elif abs(val - round(val)) < 0.001:
            row += f" {int(round(val)):>10d}"
        else:
            frac = Fraction(val).limit_denominator(10000)
            row += f" {str(frac):>10s}"
    row += f" {max_err:>10.6f}"
    print(row)

C_full = np.array(C_full)
overall_max_err = max(np.max(np.abs(Y[:, k] - X @ C_full[k])) for k in range(n))
print(f"\nOverall max residual: {overall_max_err:.6e}")
if overall_max_err < 1.0:
    print("=> EXACT linear fit confirmed (all residuals < 1)")
else:
    print(f"=> WARNING: Linear model does NOT fit exactly (max residual = {overall_max_err:.2f})")

# =====================================================================
# STEP 4: Check predicted null vectors from THM-065
# =====================================================================

print(f"\n{'=' * 72}")
print("STEP 4: Verify predicted null vectors (THM-065)")
print("=" * 72)

C_inv = C_full[:, 1:]  # (9, 7) coefficient matrix without constant
print(f"\nCoefficient matrix shape: {C_inv.shape}")

# Null vector 1 (f=4 group): (-2)*t5 + (1)*bc33
# Column order: t3=0, t5=1, t7=2, t9=3, bc33=4, bc35=5, a3=6
nv1 = np.array([0, -2, 0, 0, 1, 0, 0], dtype=float)
Cv1 = C_inv @ nv1
max_cv1 = np.max(np.abs(Cv1))

print(f"\n--- Null vector 1 (f=4 group): -2*t5 + 1*bc33 ---")
print(f"  C @ v = {[f'{x:.6f}' for x in Cv1]}")
print(f"  Max |C @ v| = {max_cv1:.2e}")
print(f"  Status: {'PASS (in null space)' if max_cv1 < 0.01 else 'FAIL'}")

combo1 = X[:, 1:] @ nv1
print(f"  Combo values: mean={np.mean(combo1):.2f}, std={np.std(combo1):.2f}, "
      f"range=[{np.min(combo1):.0f}, {np.max(combo1):.0f}]")

# Null vector 2 (f=2 group): (-2)*t7 + (1)*bc35
nv2 = np.array([0, 0, -2, 0, 0, 1, 0], dtype=float)
Cv2 = C_inv @ nv2
max_cv2 = np.max(np.abs(Cv2))

print(f"\n--- Null vector 2 (f=2 group): -2*t7 + 1*bc35 ---")
print(f"  C @ v = {[f'{x:.6f}' for x in Cv2]}")
print(f"  Max |C @ v| = {max_cv2:.2e}")
print(f"  Status: {'PASS (in null space)' if max_cv2 < 0.01 else 'FAIL'}")

combo2 = X[:, 1:] @ nv2
print(f"  Combo values: mean={np.mean(combo2):.2f}, std={np.std(combo2):.2f}, "
      f"range=[{np.min(combo2):.0f}, {np.max(combo2):.0f}]")

# Null vector 3 (f=2 group): (-4)*t7 + (1)*a3
nv3 = np.array([0, 0, -4, 0, 0, 0, 1], dtype=float)
Cv3 = C_inv @ nv3
max_cv3 = np.max(np.abs(Cv3))

print(f"\n--- Null vector 3 (f=2 group): -4*t7 + 1*a3 ---")
print(f"  C @ v = {[f'{x:.6f}' for x in Cv3]}")
print(f"  Max |C @ v| = {max_cv3:.2e}")
print(f"  Status: {'PASS (in null space)' if max_cv3 < 0.01 else 'FAIL'}")

combo3 = X[:, 1:] @ nv3
print(f"  Combo values: mean={np.mean(combo3):.2f}, std={np.std(combo3):.2f}, "
      f"range=[{np.min(combo3):.0f}, {np.max(combo3):.0f}]")

# =====================================================================
# STEP 5: Numerical rank and SVD null space analysis
# =====================================================================

print(f"\n{'=' * 72}")
print("STEP 5: SVD rank and null space dimension")
print("=" * 72)

# Use rows 0..4 (a_k = a_{n-1-k} palindromy)
C_half = C_inv[:5, :]
U, S, Vt = np.linalg.svd(C_half, full_matrices=True)

print(f"\nCoefficient matrix (rows k=0..4, 7 invariant columns):")
print(f"  Shape: {C_half.shape}")
print(f"\nSingular values:")
for i, s in enumerate(S):
    marker = " << ZERO" if s < 1e-6 * S[0] else ""
    print(f"  sigma_{i} = {s:15.6f}{marker}")

rank = np.sum(S > 1e-6 * S[0])
null_dim = C_half.shape[1] - rank

print(f"\nNumerical rank: {rank}")
print(f"Null space dimension: {null_dim}")
print(f"Expected (THM-065): rank=4, null_dim=3")
print(f"Match: {'PASS' if rank == 4 and null_dim == 3 else 'FAIL'}")

# Extract SVD null vectors
print(f"\nSVD null vectors (rows of V^T beyond rank):")
svd_null_vectors = []
for i in range(rank, Vt.shape[0]):
    v = Vt[i]
    svd_null_vectors.append(v)
    max_abs = max(abs(x) for x in v)
    v_norm = v / max_abs
    terms = []
    for j, name in enumerate(inv_names):
        if abs(v_norm[j]) > 0.01:
            frac = Fraction(v_norm[j]).limit_denominator(100)
            terms.append(f"({frac})*{name}")
    print(f"  NV{i-rank+1}: {' + '.join(terms)}")
    Cv = C_half @ v
    print(f"        Max |C @ v| = {max(abs(x) for x in Cv):.2e}")

# =====================================================================
# STEP 6: f-level weighted sums
# =====================================================================

print(f"\n{'=' * 72}")
print("STEP 6: f-level weighted sums (what a_k determines)")
print("=" * 72)

print(f"\nThe forward-edge distribution determines exactly these 4 quantities:")
print(f"  S_6 = 2^1 * t3 = 2*t3")
print(f"  S_4 = 2^1 * t5 + 2^2 * bc33 = 2*t5 + 4*bc33")
print(f"  S_2 = 2^1 * t7 + 2^2 * bc35 + 2^3 * a3 = 2*t7 + 4*bc35 + 8*a3")
print(f"  S_0 = 2^1 * t9 = 2*t9")

print(f"\nSample f-level sums:")
for trial in range(5):
    inv = X[trial, 1:]
    t3, t5, t7, t9, bc33, bc35, a3 = [int(x) for x in inv]
    S6 = 2 * t3
    S4 = 2 * t5 + 4 * bc33
    S2 = 2 * t7 + 4 * bc35 + 8 * a3
    S0 = 2 * t9
    H = 1 + S6 + S4 + S2 + S0
    print(f"  T{trial}: S6={S6:4d}, S4={S4:4d}, S2={S2:4d}, S0={S0:2d}, "
          f"H = 1+{S6}+{S4}+{S2}+{S0} = {H}")

# =====================================================================
# STEP 7: Sample tournament details
# =====================================================================

print(f"\n{'=' * 72}")
print("STEP 7: Explicit null-vector verification on sample tournaments")
print("=" * 72)

print(f"\nFor each tournament, the null combinations vary but are invisible to a_k:")
print(f"  NV1: -2*t5 + bc33       (f=4 group)")
print(f"  NV2: -2*t7 + bc35       (f=2 group)")
print(f"  NV3: -4*t7 + a3         (f=2 group)")
print()

sample_indices = list(range(0, 200, 20))
for idx in sample_indices:
    inv = X[idx, 1:]
    t3, t5, t7, t9, bc33, bc35, a3 = [int(x) for x in inv]
    nv1_val = -2 * t5 + bc33
    nv2_val = -2 * t7 + bc35
    nv3_val = -4 * t7 + a3

    pred_ak = [int(round(X[idx] @ C_full[k])) for k in range(n)]
    actual_ak = [int(Y[idx, k]) for k in range(n)]
    match = all(p == a for p, a in zip(pred_ak, actual_ak))

    print(f"  T{idx:3d}: t3={t3:2d} t5={t5:3d} t7={t7:3d} t9={t9} "
          f"bc33={bc33:3d} bc35={bc35:4d} a3={a3:2d} | "
          f"NV1={nv1_val:5d} NV2={nv2_val:5d} NV3={nv3_val:5d} | "
          f"a_k fit: {'OK' if match else 'FAIL'}")

# =====================================================================
# SUMMARY
# =====================================================================

print(f"\n{'=' * 72}")
print("SUMMARY")
print("=" * 72)

all_pass = True

if overall_max_err >= 1.0:
    print(f"[FAIL] Linear model does not fit exactly (max residual = {overall_max_err:.2f})")
    all_pass = False
else:
    print("[PASS] Exact linear fit: a_k is linear in OCF invariants")

if rank == 4 and null_dim == 3:
    print(f"[PASS] Null space dimension = {null_dim} (predicted by THM-065)")
else:
    print(f"[FAIL] Null space dimension = {null_dim}, expected 3")
    all_pass = False

for nv, name, max_cv in [(nv1, "-2*t5 + bc33", max_cv1),
                          (nv2, "-2*t7 + bc35", max_cv2),
                          (nv3, "-4*t7 + a3", max_cv3)]:
    if max_cv < 0.01:
        print(f"[PASS] Null vector {name}: max |C*v| = {max_cv:.2e}")
    else:
        print(f"[FAIL] Null vector {name}: max |C*v| = {max_cv:.2e}")
        all_pass = False

for combo, name in [(combo1, "-2*t5+bc33"), (combo2, "-2*t7+bc35"), (combo3, "-4*t7+a3")]:
    if np.std(combo) > 1.0:
        print(f"[PASS] {name} is nontrivial (std={np.std(combo):.2f})")
    else:
        print(f"[WARN] {name} has low variance (std={np.std(combo):.2f})")

if ocf_errors > 0:
    print(f"[FAIL] OCF failed for {ocf_errors}/{NUM_TOURNAMENTS} tournaments")
    all_pass = False
else:
    print(f"[PASS] OCF verified for all {NUM_TOURNAMENTS} tournaments")

print(f"\nOverall: {'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
print(f"THM-065 predictions at n=9: {'CONFIRMED' if all_pass else 'NOT CONFIRMED'}")

print(f"\n{'=' * 72}")
print("DONE")
print("=" * 72)
